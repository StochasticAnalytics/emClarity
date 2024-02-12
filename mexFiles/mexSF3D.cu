
#include "include/core_headers.cuh"



#define MAX_EPSILON_ERROR 5e-3f



////////////////////////////////////////////////////////////////////////////////
// Constants


// Texture reference for 2D float texture
// texture<float, 2, cudaReadModeElementType> tex;
const float  slice_thickness_pixel_radius = 3; // TODO should
const float cosine_edge_arg = PI / (float)slice_thickness_pixel_radius;
const float cosine_edge_norm = 1.0f;// / 3.0f;


////////////////////////////////////////////////////////////////////////////////
//! Transform an image using texture lookups
//! @param outputData  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__global__ void sf3dKernel(const cudaTextureObject_t tex_obj,
                           float *outputData,
                           float3 size_shift,
                           uint3 dims,
                           float2 sinAcosA
                           ) {

  // Assuming a single-Y-axis tilt such that the rotation is
  // [ c, 0, s,
  //   0, 1, 0,
  //  -s, 0, c]

  int x = blockIdx.x*blockDim.x + threadIdx.x;
  if (x >= dims.x) { return ; }
  int y = blockIdx.y*blockDim.y + threadIdx.y;
  if (y >= dims.y) { return ; }

  // For every 3D voxel invert the tilt transformation to see where the data would be in the xy plane of the 2D image
  float tu,tw,tw_pre,zWeight;


  // The centered coordinate in the 3d volume
  float x_centered = (float)x - (float)dims.x/2;
  tw_pre = -x_centered * sinAcosA.x + size_shift.z;
  float y_centered_normalized = ((float)y - (float)(dims.y/2) + size_shift.y) / (float)dims.y + 0.5f;

  for (int z = 0; z < dims.z; z++) {

    float z_centered = float(z - dims.z/2);
    tw = tw_pre + z_centered * sinAcosA.y;

    if (tw < -slice_thickness_pixel_radius || tw > slice_thickness_pixel_radius) {
      continue;
    }

    // FIXME this should approximate a sinc
    zWeight = (0.5 + 0.5*cosf(tw * cosine_edge_arg )) * cosine_edge_norm;

    tu =  x_centered*sinAcosA.y + z_centered*sinAcosA.x + size_shift.x; 
    tu /= (float)dims.x; // Normalized coords
    tw /= (float)dims.z;
      

       
      outputData[ (z*dims.y + y) * dims.x + x ] += ( zWeight * tex2D<float>(tex_obj, tu + 0.5f, y_centered_normalized));

      // outputData[ (z*dims.y + y) * dims.x + x ] += ( zWeight * tex2D(tex, tu, tv));
  }

}


////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{


//  // Initialize the MathWorks GPU API.
//  mxInitGPU();

  // TODO if these aren't typed properly in the mex call, I can't cast them appropriatley here
  // TODO general angles, here assume single-Y-axis tilt
  /* Check for proper number of arguments. TODO add checks on narg and types*/
    if ( nrhs != 15 & nrhs != 16 ) 
    {
        mexErrMsgIdAndTxt("MATLAB:mexSF3D:rhs",
            "requires 14 inputs.\n");
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
  float* wiener_constant = (float *) mxGetData(prhs[15]);

  float * d_output_img = NULL;
  float * d_ctf_img = NULL;
  uint3 dims;
  uint2 ctf_dims;
  uint2 o_ctf_dims;

  float2 sinAcosA[*nTilts];



  mxGPUArray * outputArray;

  float sin_t = 0.0f;
  float cos_t = 0.0f;
  for (int iAng = 0; iAng < *nTilts; iAng++) 
  {

    sincosf( deg_2_rad(-tiltAngle[iAng]), &sin_t, &cos_t);
    sinAcosA[iAng] = make_float2(sin_t, cos_t);

  }
    

  dims     = make_uint3(wantedSize[0],wantedSize[1],wantedSize[2]);
  ctf_dims = make_uint2(wantedSize[0],wantedSize[1]);



  if (*doHalfGrid )
  {
    o_ctf_dims = make_uint2(0, ctf_dims.y/2);
    ctf_dims.x = ctf_dims.x/2 + 1;
  }
  else
  {
    o_ctf_dims = make_uint2(ctf_dims.x/2, ctf_dims.y/2);
  }




  // Allocate device memory for result
  mwSize output_dims = 3;
  mwSize output_size[3] = {dims.x, dims.y, dims.z};

  // Allocate device memory for the weights
  // Create MX array and init with zeros
  outputArray = mxGPUCreateGPUArray(output_dims,
                                    output_size,
                                    mxSINGLE_CLASS,
                                    mxREAL,
                                    MX_GPU_INITIALIZE_VALUES);

  d_output_img = (float *)(mxGPUGetData(outputArray));
  cudaTextureObject_t tex_obj = 0;
  cudaArray* cuArray = 0;

  ////////////////////////////////////////////////////////////////////


  // Allocate array and copy image data
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  checkCudaErrors(cudaMallocArray(&cuArray, &channelDesc, ctf_dims.x, ctf_dims.y));


  // TODO checkout cudaCreateChannelDescHalf https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#sixteen-bit-floating-point-textures
  struct cudaResourceDesc resDesc;
  (memset(&resDesc, 0, sizeof(resDesc)));
  resDesc.resType         = cudaResourceTypeArray;
  resDesc.res.array.array = cuArray;

  struct cudaTextureDesc texDesc;
  (memset(&texDesc, 0, sizeof(texDesc)));

  texDesc.filterMode       = cudaFilterModeLinear;
  texDesc.readMode         = cudaReadModeElementType;
  texDesc.normalizedCoords = 1;
  texDesc.addressMode[0]   = cudaAddressModeBorder; 
  texDesc.addressMode[1]   = cudaAddressModeBorder; 
  texDesc.addressMode[2]   = cudaAddressModeBorder; 




  size_t n_to_copy = ctf_dims.x * ctf_dims.y * sizeof(float);
  checkCudaErrors(cudaMallocAsync(&d_ctf_img, n_to_copy, cudaStreamPerThread));
  // d_ctf_img = (float *)(mxGPUGetData(ctfArray));



  // Would specifying a 3d grid speed up by improving drop out over a block?
  
  dim3 threads_per_block = dim3(32,32,1); // max is 1024 threads/block for 2.x --> 7.5 compute capability
  dim3 dimGrid = dim3((dims.x + threads_per_block.x -1) / threads_per_block.x,
                      (dims.y + threads_per_block.y -1) / threads_per_block.y,
                      1);//(dims.z + threads_per_block.z -1) / threads_per_block.z);

  dim3 ctfBlock(32, 32, 1);
  dim3 ctfGrid((ctf_dims.x + ctfBlock.x - 1) / ctfBlock.x, (ctf_dims.y + ctfBlock.y - 1) / ctfBlock.y, 1);
  bool calc_centered = true;

  float3 size_shift = make_float3(0.0f, 0.0f, 0.0f);
  if (IsEven(dims.x)) size_shift.x = 0.5f;
  if (IsEven(dims.y)) size_shift.y = 0.5f;
  if (IsEven(dims.z)) size_shift.z = 0.5f;

  // Calculate this prior to any half-dim reduction
  float2 fourierVoxelSize;
  fourierVoxelSize = make_float2( 1.0f/(pixelSize[0] * (float)ctf_dims.x), 
                                  1.0f/(pixelSize[0] * (float)ctf_dims.y));

  bool is_initialized = false;

  for (int iAng = 0 ; iAng < *nTilts ; iAng ++)
  {

    if (iAng > 0) {
      fourierVoxelSize.x = 1.0f/(pixelSize[iAng] * (float)ctf_dims.x);
      fourierVoxelSize.y = 1.0f/(pixelSize[iAng] * (float)ctf_dims.y);
    }


    ctfParams b_ctf(*doHalfGrid,*doSqCTF,pixelSize[iAng],waveLength[iAng],CS[iAng],AmpContrast[iAng],
                    defocus1[iAng],  defocus2[iAng], defocusAst[iAng]);


    // Create the 2d ctf
    ctf<<< ctfGrid, ctfBlock ,0,cudaStreamPerThread >>>(d_ctf_img, ctf_dims, o_ctf_dims, b_ctf, fourierVoxelSize,
                                 calc_centered, occupancy[iAng], exposure[iAng], *wiener_constant);


    // FIXME
    checkCudaErrors(cudaMemcpyToArrayAsync(cuArray, 0, 0, d_ctf_img, n_to_copy, cudaMemcpyDeviceToDevice, cudaStreamPerThread));
    // checkCudaErrors(cudaMemcpy2DToArrayAsync(cuArray, 0, 0, d_ctf_img, ctf_dims.x * sizeof(float), ctf_dims.x, ctf_dims.y, cudaMemcpyDeviceToDevice, cudaStreamPerThread));

    if (! is_initialized) {
        cudaStreamSynchronize(cudaStreamPerThread);
        checkCudaErrors(cudaCreateTextureObject(&tex_obj, &resDesc, &texDesc, NULL));
        is_initialized = true;  
    }
    // Call the sf3d kernel
    sf3dKernel<<<dimGrid, threads_per_block,0,cudaStreamPerThread >>>(tex_obj, d_output_img, size_shift, dims, sinAcosA[iAng]);

  }

  checkCudaErrors(cudaStreamSynchronize(cudaStreamPerThread));

  // Check if kernel execution generated an error
  getLastCudaError("Kernel execution failed\n");


  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);


  checkCudaErrors(cudaFreeArray(cuArray));
  checkCudaErrors(cudaDestroyTextureObject(tex_obj));

  mxGPUDestroyGPUArray(outputArray);

  // mxGPUDestroyGPUArray(ctfArray);
  cudaFreeAsync(d_ctf_img, cudaStreamPerThread);


//  (cudaFree(d_input_img));
//  (cudaFree(d_output_img));


}
