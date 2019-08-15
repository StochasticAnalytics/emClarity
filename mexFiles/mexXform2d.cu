#include "include/core_headers.cuh"



#define MAX_EPSILON_ERROR 5e-3f



////////////////////////////////////////////////////////////////////////////////
// Constants


// Texture reference for 2D float texture
texture<float, 2, cudaReadModeElementType> tex;
const float EXTRAPVAL = 0.0f;


////////////////////////////////////////////////////////////////////////////////
//! Transform an image using texture lookups
//! @param outputData  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__global__ void transformKernel_FWD(float *outputData,
                                uint2 dims,
                                float2 rm_1,
                                float2 rm_2,
                                float2 shifts,
                                float extrapVal,
                                bool doFwdXform)
{


  // calculate normalized texture coordinates
  unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
  float u,v,tu,tv;

  if (doFwdXform)
  {
    u = (float)x - (float)dims.x/2; 
    v = (float)y - (float)dims.y/2; 
    tu = u*rm_1.x + v*rm_1.y + shifts.x; 
    tv = u*rm_2.x + v*rm_2.y + shifts.y;
  }
  else  
  {
    u = (float)x - (float)dims.x/2 - shifts.x; 
    v = (float)y - (float)dims.y/2 - shifts.y; 
    tu = u*rm_1.x + v*rm_1.y; 
    tv = u*rm_2.x + v*rm_2.y;
  }


  tu /= (float)dims.x; 
  tv /= (float)dims.y; 
  tu += 0.5f;
  tv += 0.5f;

  // TODO one of these is probably supposed to be inclusive
  if (tu < 0 | tv < 0 | tu > 1 | tv > 1)
  {
    outputData[y*dims.x + x] = extrapVal;
  }
  else
  {
    outputData[y*dims.x + x] = tex2D(tex, tu, tv);
  }

}


////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{

  int iSpot = 1;
  float * d_output_img = NULL;
  float * d_input_img  = NULL;
  uint2 dims;
  float2 shifts;
  float2 rm_1;
  float2 rm_2;
  bool *doFwdXform; // rotate then shift
  unsigned int size;

/* Check for proper number of arguments. TODO add checks on narg and types*/
  if ( ! mxGPUIsValidGPUData(prhs[0]) ) 
  {
      mexErrMsgIdAndTxt("MATLAB:mexFFT:rhs",
          "This inputArray is not valid gpudata.");
  }

  

  mxGPUArray const * inputArray  = mxGPUCreateFromMxArray(prhs[0]);
  mxGPUArray * outputArray;

  float *angles = (float *) mxGetData(prhs[1]);
  float *ts = (float *) mxGetData(prhs[2]);
  doFwdXform = (bool *) mxGetData(prhs[3]);

  if (*doFwdXform)
  {
    mexPrintf("doing forward\n");
    // transpose matrix
    shifts = make_float2(-ts[0],-ts[1]);
    rm_1   = make_float2(angles[0],angles[1]);
    rm_2   = make_float2(angles[3],angles[4]);
    mexPrintf("%f %f %f %f shifts %f %f\n",rm_1.x,rm_1.y,rm_2.x,rm_2.y, shifts.x, shifts.y);
    
  }
  else
  {
    mexPrintf("doing inv\n");

    shifts = make_float2(ts[0],ts[1]);
    rm_1   = make_float2(angles[0],angles[3]);
    rm_2   = make_float2(angles[1],angles[4]);
  }

  d_input_img = (float *)(mxGPUGetDataReadOnly(inputArray));

  mwSize const   input_dims = mxGPUGetNumberOfDimensions(inputArray);
  mwSize const * input_size = mxGPUGetDimensions(inputArray);
  mwSize const   numel_input = mxGPUGetNumberOfElements(inputArray);
  mxComplexity   input_data_type = mxGPUGetComplexity(inputArray);

  dims = make_uint2(input_size[0],input_size[1]);
  

  size = dims.x * dims.y * sizeof(float);
  // Allocate device memory for result
  float *dData = NULL;

  // Create MX array and init with zeros
  outputArray = mxGPUCreateGPUArray(input_dims,
                                    input_size,
                                    mxSINGLE_CLASS,
                                    input_data_type,
                                    MX_GPU_INITIALIZE_VALUES);

  d_output_img = (float *)(mxGPUGetData(outputArray));

  cudaArray *cuArray;
  // TODO where does the 32 come from?
  cudaChannelFormatDesc channelDesc =
        cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

  checkCudaErrors(cudaMallocArray(&cuArray,
                                  &channelDesc,
                                  dims.x,
                                  dims.y));
  checkCudaErrors(cudaMemcpyToArray(cuArray,
                                    0,
                                    0,
                                    d_input_img,
                                    size,
                                    cudaMemcpyDeviceToDevice));

  // Set texture parameters
  // cudaAddressModeWrap cudaAddressModeClamp cudaAddressModeMirror cudaAddressModeBorder
  tex.addressMode[0] = cudaAddressModeClamp;
  tex.addressMode[1] = cudaAddressModeClamp;

  tex.filterMode = cudaFilterModeLinear;

  tex.normalized = true;    // access with normalized texture coordinates



  
  // Bind the array to the texture
  checkCudaErrors(cudaBindTextureToArray(tex, cuArray, channelDesc));


  // How to choose these?
  dim3 dimBlock(8, 8, 1);

  dim3 dimGrid(dims.x / dimBlock.x, dims.y / dimBlock.y, 1);

  


  transformKernel_FWD<<<dimGrid, dimBlock, 0>>>(d_output_img, dims,
                                              rm_1,rm_2,shifts,EXTRAPVAL,*doFwdXform);
  




  checkCudaErrors(cudaDeviceSynchronize());


  // Check if kernel execution generated an error
  getLastCudaError("Kernel execution failed");


  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);

  mxGPUDestroyGPUArray(inputArray);
  mxGPUDestroyGPUArray(outputArray);


//  checkCudaErrors(cudaFree(d_input_img));
//  checkCudaErrors(cudaFree(d_output_img));


}
