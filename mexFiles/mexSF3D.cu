
#include "include/core_headers.cuh"



#define MAX_EPSILON_ERROR 5e-3f



////////////////////////////////////////////////////////////////////////////////
// Constants


// Texture reference for 2D float texture
texture<float, 2, cudaReadModeElementType> tex;
const float EXTRAPVAL = 0.0f;
const int  slice_thickness_pixel_radius = 5; // TODO should

const float wanted_padding = 1.0; // oversample the 2d ctf

////////////////////////////////////////////////////////////////////////////////
//! Transform an image using texture lookups
//! @param outputData  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__global__ void sf3dKernel(float *outputData,
                           uint3 dims,
                           float2 sinAcosA,
                           float extrapVal)
{

  // Assuming a single-Y-axis tilt such that the rotation is
  // [ c, 0, s,
  //   0, 1, 0,
  //  -s, 0, c]

  unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
  if (x >= dims.x) { return ; }
  unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
  if (y >= dims.y) { return ; }
  unsigned int z = blockIdx.z; //*blockDim.z + threadIdx.z;
  if (z >= dims.z) {  return ; }
  float u,w,tv,tu,tw,zWeight;


  // First calc the Z-dimension, and check that we are close enough to the plane
  u = (float)x - (float)dims.x/2;
  w = (float)z - (float)dims.z/2;
  tw = -u*sinAcosA.x + w*sinAcosA.y;

  if (tw >= -slice_thickness_pixel_radius & tw <= slice_thickness_pixel_radius)
  {
    // FIXME this should approximate a sinc
    zWeight = 0.5 + 0.5*cosf(tw*PI/(float)(slice_thickness_pixel_radius+1.0f));

    tu =  u*sinAcosA.y + w*sinAcosA.x; 
    tu /= (float)dims.x; // Normalized coords
    tw /= (float)dims.z;

    // re-use u to calc a radial weight. Set u at origin to u(1) as in imod
    if (u == 0) { u = 0.2; }
    
    // FIXME this extra division could go
//    u /= (float)(dims.x/2);
    u /= 100.0f;

    tu += 0.5f;
    tw += 0.5f;
    tv = ((float)y - (float)dims.y/2) / (float)dims.y + 0.5f; 

    // TODO one of these is probably supposed to be inclusive
    if (tu > 0 & tw > 0 & tu > 0 & tw < 1 & tv < 1 & tv < 1)
    {
      // TODO The radial weighting and exposure weighting can, and probably should just be done on the 2d ctf prior to texturing
      outputData[ (z*dims.y + y) * dims.x + x ] += ( zWeight * (fabsf(u)) * tex2D(tex, tu, tv));
    }
  }

}


////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{


  // Initialize the MathWorks GPU API.
  mxInitGPU();

  // TODO if these aren't typed properly in the mex call, I can't cast them appropriatley here
  // TODO general angles, here assume single-Y-axis tilt
  /* Check for proper number of arguments. TODO add checks on narg and types*/
    if ( nrhs != 14 & nrhs != 15 ) 
    {
        mexErrMsgIdAndTxt("MATLAB:mexSF3D:rhs",
            "requires 14 inputs.");
    }

  bool * doHalfGrid =  (bool *) mxGetData(prhs[0]);  
  bool * doSqCTF    =  (bool *) mxGetData(prhs[1]); 
  uint * wantedSize =  (uint *) mxGetData(prhs[2]);
  float* pixelSize  =  (float*) mxGetData(prhs[3]);
  float* waveLength =  (float*) mxGetData(prhs[4]);
  float* CS         =  (float*) mxGetData(prhs[5]);
  float* defocus1   =  (float*) mxGetData(prhs[6]);
  float* defocus2   =  (float*) mxGetData(prhs[7]);
  float* defocusAst =  (float*) mxGetData(prhs[8]);
  float* AmpContrast=  (float*) mxGetData(prhs[9]);
  uint*   nTilts     = (uint  *) mxGetData(prhs[10]);
  float* tiltAngle  =  (float*) mxGetData(prhs[11]);
  float* exposure   =  (float*) mxGetData(prhs[12]);
  float* occupancy  =  (float*) mxGetData(prhs[13]);
  int*   launch     = (int *) mxGetData(prhs[14]);

  float * d_output_img = NULL;
  float * d_ctf_img = NULL;
  uint3 dims;
  uint2 ctf_dims;
  uint2 o_ctf_dims;

  float2 sinAcosA[*nTilts];


  ctfParams b_ctf(*doHalfGrid,*doSqCTF,*pixelSize,*waveLength,*CS,*AmpContrast,
                  *defocus1,  *defocus2, *defocusAst);

///* Check for proper number of arguments. TODO add checks on narg and types*/
//  if ( ! mxGPUIsValidGPUData(prhs[0]) ) 
//  {
//      mexErrMsgIdAndTxt("MATLAB:mexFFT:rhs",
//          "This inputArray is not valid gpudata.");
//  }


  mxGPUArray * outputArray;

  for (int iAng = 0; iAng < *nTilts; iAng++) 
  {
    sinAcosA[iAng] = make_float2(sinf(tiltAngle[iAng]*PI/180.0f),cosf(tiltAngle[iAng]*PI/180.0f));

  }
    

  dims     = make_uint3(wantedSize[0],wantedSize[1],wantedSize[2]);
  ctf_dims = make_uint2(wantedSize[0]*wanted_padding,wantedSize[1]*wanted_padding);

  // Calculate this prior to any half-dim reduction
  float2 fourierVoxelSize;
  fourierVoxelSize = make_float2( 1.0f/(*pixelSize * (float)ctf_dims.x), 
                                  1.0f/(*pixelSize * (float)ctf_dims.y));

  if (*doHalfGrid )
  {
    o_ctf_dims = make_uint2(0, ctf_dims.y/2);
    ctf_dims.x = ctf_dims.x/2 + 1;
  }
  else
  {
    o_ctf_dims = make_uint2(ctf_dims.x/2, ctf_dims.y/2);
  }



  long numel_output, numel_ctf;
  numel_output = dims.x * dims.y  *dims.z * sizeof(float);
  numel_ctf    = ctf_dims.x * ctf_dims.y * sizeof(float);

  // Allocate device memory for result
  mwSize output_dims = 3;
  mwSize output_size[3] = {dims.x, dims.y, dims.z};

  // Create MX array and init with zeros
  outputArray = mxGPUCreateGPUArray(output_dims,
                                    output_size,
                                    mxSINGLE_CLASS,
                                    mxREAL,
                                    MX_GPU_INITIALIZE_VALUES);

  d_output_img = (float *)(mxGPUGetData(outputArray));

  cudaArray *cuArray;

  // TODO where does the 32 come from?
  cudaChannelFormatDesc channelDesc =
        cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);


//  mexit(__LINE__);
  checkCudaErrors(cudaMallocArray(&cuArray,
                                  &channelDesc,
                                  ctf_dims.x,
                                  ctf_dims.y));

  // Set texture parameters
  // cudaAddressModeWrap cudaAddressModeClamp cudaAddressModeMirror cudaAddressModeBorder
  tex.addressMode[0] = cudaAddressModeClamp;
  tex.addressMode[1] = cudaAddressModeClamp;

  tex.filterMode = cudaFilterModeLinear;

  tex.normalized = true;    // access with normalized texture coordinates

  // Params for the 2d ctf
  mwSize const   ctf_number_of_dims = 2;
  mxComplexity   ctf_type = mxREAL;
  mxClassID      ctf_class = mxSINGLE_CLASS;
  mwSize         ctf_size[2];
  ctf_size[0] = (mwSize)  ctf_dims.x;
  ctf_size[1] = (mwSize)  ctf_dims.y;
    
  mxGPUArray *ctfArray;  
  // TODO it would be nice not to init all the zeros, but then the fourier padding would need to be dealt with.
  ctfArray = mxGPUCreateGPUArray(ctf_number_of_dims,
                                 ctf_size,
                                 ctf_class,
                                 ctf_type,
                                 MX_GPU_INITIALIZE_VALUES);

  d_ctf_img = (cufftReal *)(mxGPUGetData(outputArray));

  dim3 ctfBlock(32, 32, 1);
  dim3 ctfGrid(ctf_dims.x / ctfBlock.x, ctf_dims.y / ctfBlock.y, 1);
  bool calc_centered = true;
  for (int iAng = 0 ; iAng < *nTilts ; iAng ++)
  {

    // Create the 2d ctf
    ctf<<< ctfGrid, ctfBlock >>>(d_ctf_img, ctf_dims, o_ctf_dims, b_ctf, fourierVoxelSize, calc_centered );


    // Put the ctf in tex2
    checkCudaErrors(cudaMemcpyToArray(cuArray,
                                      0,
                                      0,
                                      d_ctf_img,
                                      numel_ctf,
                                      cudaMemcpyDeviceToDevice));

    
    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(tex, cuArray, channelDesc));

    // Call the sf3d kernel
    //TODO

  // Would specifying a 3d grid speed up by improving drop out over a block?
    int dimDist = 32;
    dim3 threads_per_block = dim3(dimDist,dimDist,1); // max is 1024 threads/block for 2.x --> 7.5 compute capability
    dim3 dimGrid = dim3((dims.x+dimDist-1) / dimDist,(dims.y+dimDist-1)/dimDist,dims.z);

    sf3dKernel<<<dimGrid, threads_per_block >>>(d_output_img, dims,
                                                  sinAcosA[iAng],EXTRAPVAL);
  }




  checkCudaErrors(cudaDeviceSynchronize());


  // Check if kernel execution generated an error
  getLastCudaError("Kernel execution failed");


  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);


  mxGPUDestroyGPUArray(outputArray);


//  checkCudaErrors(cudaFree(d_input_img));
//  checkCudaErrors(cudaFree(d_output_img));


}
