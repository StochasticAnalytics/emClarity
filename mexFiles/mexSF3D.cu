
#include "include/core_headers.cuh"



#define MAX_EPSILON_ERROR 5e-3f



////////////////////////////////////////////////////////////////////////////////
// Constants


// Texture reference for 2D float texture
texture<float, 2, cudaReadModeElementType> tex;
const float EXTRAPVAL = 0.0f;
const float  slice_thickness_pixel_radius = 4; // TODO should
const float cosine_edge_arg = PI / (float)slice_thickness_pixel_radius;
const float cosine_edge_norm = 1.0f / 4.0f;

const float wanted_padding = 1.0; // oversample the 2d ctf

////////////////////////////////////////////////////////////////////////////////
//! Transform an image using texture lookups
//! @param outputData  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__global__ void sf3dKernel(cudaTextureObject_t* thisTexObj,
                           uint nTilts,
                           float *outputData,
                           float3 size_shift,
                           uint3 dims,
                           float2* sinAcosA,
                           float extrapVal)
{

  // Assuming a single-Y-axis tilt such that the rotation is
  // [ c, 0, s,
  //   0, 1, 0,
  //  -s, 0, c]

  int x = blockIdx.x*blockDim.x + threadIdx.x;
  if (x >= dims.x) { return ; }
  int y = blockIdx.y*blockDim.y + threadIdx.y;
  if (y >= dims.y) { return ; }
  int z = blockIdx.z*blockDim.z + threadIdx.z;
  if (z >= dims.z) {  return ; }
  float u,w,tv,tu,tw,zWeight;


  // First calc the Z-dimension, and check that we are close enough to the plane
  u = (float)x - (float)dims.x/2;
  w = (float)z - (float)dims.z/2;

  tv = ((float)y - (float)dims.y/2 + size_shift.y) / (float)dims.y + 0.5f; 
  int idx = (z*dims.y + y) * dims.x + x;

  for (int iAng = 0 ; iAng < (int)nTilts ; iAng ++)
  {
    tw = -u*sinAcosA[iAng].x + w*sinAcosA[iAng].y + size_shift.z;

    if (tw >= -slice_thickness_pixel_radius & tw <= slice_thickness_pixel_radius)
    {
      // FIXME this should approximate a sinc
      zWeight = (0.5 + 0.5*cosf(tw * cosine_edge_arg )) * cosine_edge_norm;

      tu =  u*sinAcosA[iAng].y + w*sinAcosA[iAng].x + size_shift.x; 
      tu /= (float)dims.x; // Normalized coords
      tw /= (float)dims.z;
      
      // FIXME this extra division could go
  //    u /= (float)(dims.x/2);
      //100.0f;

      tu += 0.5f;
      tw += 0.5f;


      // TODO one of these is probably supposed to be inclusive
      if (tu > 0 & tw > 0 & tv > 0 & tu < 1 - 1/(float)dims.x & tw < 1 - 1/(float)dims.z)
      {
  //      // re-use u to calc a radial weight. Set u at origin to u(1) as in imod
  //      if (u == 0) { u = 0.2; }
  //      u /= (float)dims.x ;


        // TODO The radial weighting and exposure weighting can, and probably should just be done on the 2d ctf prior to texturing
  //      outputData[ idx ] += ( zWeight * (fabsf(u)) * tex2D(tex, tu, tv));
        outputData[ idx ] += ( zWeight * tex2D<float>(thisTexObj[iAng], tu, tv) );
      }
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
  uint * wantedSize =  (uint *) mxGetData(prhs[2]); // Should be a uint32 in matlab
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
  int*   launch     = (int *) mxGetData(prhs[14]); // should be an int16 in matlab

  float * d_output_img = NULL;
  float * d_ctf_img = NULL;
  uint3 dims;
  uint2 ctf_dims;
  uint2 o_ctf_dims;

  float2 sinAcosA[*nTilts];

//	cudaStream_t calcStream;
//	cudaEvent_t  calcEvent;
  cudaTextureObject_t tex[*nTilts];
  cudaArray_t cuArray[*nTilts];

  checkCudaErrors(cudaMalloc((void **)&tex, *nTilts * sizeof(cudaTextureObject_t)));
  checkCudaErrors(cudaMalloc((void **)&cuArray, *nTilts * sizeof(cudaTextureObject_t)));




///* Check for proper number of arguments. TODO add checks on narg and types*/
//  if ( ! mxGPUIsValidGPUData(prhs[0]) ) 
//  {
//      mexErrMsgIdAndTxt("MATLAB:mexFFT:rhs",
//          "This inputArray is not valid gpudata.");
//  }

mexPrintf("%d\n",__LINE__);
  mxGPUArray * outputArray;
mexPrintf("%d\n",__LINE__);

  for (int iAng = 0; iAng < *nTilts; iAng++) 
  {

    sinAcosA[iAng] = make_float2(sinf(-tiltAngle[iAng]*PI/180.0f),cosf(-tiltAngle[iAng]*PI/180.0f));

  }
    

  dims     = make_uint3(wantedSize[0],wantedSize[1],wantedSize[2]);
  ctf_dims = make_uint2(wantedSize[0]*wanted_padding,wantedSize[1]*wanted_padding);

mexPrintf("%d\n",__LINE__);

  if (*doHalfGrid )
  {
    o_ctf_dims = make_uint2(0, ctf_dims.y/2);
    ctf_dims.x = ctf_dims.x/2 + 1;
  }
  else
  {
    o_ctf_dims = make_uint2(ctf_dims.x/2, ctf_dims.y/2);
  }


mexPrintf("%d\n",__LINE__);
  long  numel_ctf;
  numel_ctf    = ctf_dims.x * ctf_dims.y * sizeof(float);
mexPrintf("%d\n",__LINE__);
  // Allocate device memory for result
  mwSize output_dims = 3;
  mwSize output_size[3] = {dims.x, dims.y, dims.z};

  // Allocate device memory for the weights
mexPrintf("%d\n",__LINE__);
  // Create MX array and init with zeros
  outputArray = mxGPUCreateGPUArray(output_dims,
                                    output_size,
                                    mxSINGLE_CLASS,
                                    mxREAL,
                                    MX_GPU_INITIALIZE_VALUES);
mexPrintf("%d\n",__LINE__);
  d_output_img = (float *)(mxGPUGetData(outputArray));



//  // TODO where does the 32 come from?
//  cudaChannelFormatDesc channelDesc =
//        cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

mexPrintf("%d\n",__LINE__);
////  mexit(__LINE__);
//  checkCudaErrors(cudaMallocArray(&cuArray,
//                                  &channelDesc,
//                                  ctf_dims.x,
//                                  ctf_dims.y));

//  // Set texture parameters
//  // cudaAddressModeWrap cudaAddressModeClamp cudaAddressModeMirror cudaAddressModeBorder
//  tex.addressMode[0] = cudaAddressModeClamp;
//  tex.addressMode[1] = cudaAddressModeClamp;

//  tex.filterMode = cudaFilterModeLinear;

//  tex.normalized = true;    // access with normalized texture coordinates

  // Params for the 2d ctf
  mwSize const   ctf_number_of_dims = 2;
  mxComplexity   ctf_type = mxREAL;
  mxClassID      ctf_class = mxSINGLE_CLASS;
  mwSize         ctf_size[2];
  ctf_size[0] = (mwSize)  ctf_dims.x;
  ctf_size[1] = (mwSize)  ctf_dims.y;
 mexPrintf("%d\n",__LINE__);   
  mxGPUArray *ctfArray;  
  // TODO it would be nice not to init all the zeros, but then the fourier padding would need to be dealt with.
  ctfArray = mxGPUCreateGPUArray(ctf_number_of_dims,
                                 ctf_size,
                                 ctf_class,
                                 ctf_type,
                                 MX_GPU_INITIALIZE_VALUES);

  d_ctf_img = (float *)(mxGPUGetData(ctfArray));



mexPrintf("%d\n",__LINE__);


  // Would specifying a 3d grid speed up by improving drop out over a block?
    int dimDist = 32;
    dim3 threads_per_block = dim3(dimDist,dimDist,1); // max is 1024 threads/block for 2.x --> 7.5 compute capability
    dim3 dimGrid = dim3((dims.x+dimDist-1) / dimDist,(dims.y+dimDist-1)/dimDist,dims.z);

    dim3 ctfBlock(32, 32, 1);
    dim3 ctfGrid((ctf_dims.x + ctfBlock.x - 1) / ctfBlock.x, (ctf_dims.y + ctfBlock.y - 1) / ctfBlock.y, 1);
    bool calc_centered = true;

  float3 size_shift = make_float3(0.0f, 0.0f, 0.0f);
  if (IsEven(dims.x)) size_shift.x = 0.5f;
  if (IsEven(dims.y)) size_shift.y = 0.5f;
  if (IsEven(dims.z)) size_shift.z = 0.5f;



 mexPrintf("%d\n",__LINE__);


    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    
mexPrintf("%d\n",__LINE__);
    cudaMemcpy3DParms p = {0};
    p.extent = make_cudaExtent(ctf_dims.x,ctf_dims.y,1);
    p.srcPtr = make_cudaPitchedPtr(d_ctf_img, ctf_dims.x*sizeof(float),ctf_dims.x,ctf_dims.y);
    p.kind = cudaMemcpyDeviceToDevice;
mexPrintf("%d\n",__LINE__);

    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
 mexPrintf("%d\n",__LINE__);   
   
    struct cudaTextureDesc texDesc;
    memset(&texDesc,0,sizeof(texDesc));
mexPrintf("%d\n",__LINE__);
    texDesc.filterMode = cudaFilterModeLinear;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = 1;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.addressMode[1] = cudaAddressModeClamp;
    texDesc.addressMode[2] = cudaAddressModeClamp;
mexPrintf("%d\n",__LINE__);
    

  for (int iAng = 0 ; iAng < *nTilts ; iAng ++)
  {

    // Calculate this prior to any half-dim reduction
    float2 fourierVoxelSize;
    fourierVoxelSize = make_float2( 1.0f/(pixelSize[iAng] * (float)ctf_dims.x), 
                                    1.0f/(pixelSize[iAng] * (float)ctf_dims.y));

mexPrintf("%d\n",__LINE__);

    ctfParams b_ctf(*doHalfGrid,*doSqCTF,pixelSize[iAng],waveLength[iAng],CS[iAng],AmpContrast[iAng],
                    defocus1[iAng],  defocus2[iAng], defocusAst[iAng]);

mexPrintf("%d\n",__LINE__);
    // Create the 2d ctf
    ctf<<< ctfGrid, ctfBlock, 0, cudaStreamPerThread >>>(d_ctf_img, ctf_dims, o_ctf_dims, b_ctf, fourierVoxelSize,
                                 calc_centered, occupancy[iAng], exposure[iAng]);


mexPrintf("%d\n",__LINE__);
    tex[iAng] = 0;
    checkCudaErrors(cudaMalloc3DArray(&cuArray[iAng],
                                      &channelDesc,
                                      make_cudaExtent(dims.x,dims.y,0)));
mexPrintf("%d\n",__LINE__);
    p.dstArray = cuArray[iAng];
    resDesc.res.array.array = cuArray[iAng];
    cudaCreateTextureObject(&tex[iAng],&resDesc,&texDesc,NULL);
    cudaMemcpy3DAsync(&p, cudaStreamPerThread);


mexPrintf("%d\n",__LINE__);

    // FIXME if you could bin an array of texture objects, you could launch the sf3dKernel outside the loop once. 
  }
mexPrintf("%d\n",__LINE__);

    // Call the sf3d kernel
    sf3dKernel<<<dimGrid, threads_per_block,0,cudaStreamPerThread>>>(tex, *nTilts, d_output_img, size_shift,
                                                                     dims, sinAcosA,EXTRAPVAL);

  cudaStreamSynchronize(cudaStreamPerThread);


mexPrintf("%d\n",__LINE__);
  // Check if kernel execution generated an error
  getLastCudaError("Kernel execution failed");

mexPrintf("%d\n",__LINE__);
  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);
mexPrintf("%d\n",__LINE__);
  for (int iAng = 0 ; iAng < *nTilts ; iAng ++)
  {
    checkCudaErrors(cudaFreeArray(cuArray[iAng]));
    checkCudaErrors(cudaDestroyTextureObject(tex[iAng]));
  }
mexPrintf("%d\n",__LINE__);
  mxGPUDestroyGPUArray(outputArray);
mexPrintf("%d\n",__LINE__);
  mxGPUDestroyGPUArray(ctfArray);

mexPrintf("%d\n",__LINE__);
//  checkCudaErrors(cudaFree(d_input_img));
//  checkCudaErrors(cudaFree(d_output_img));


}
