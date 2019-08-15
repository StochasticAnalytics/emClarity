#include "include/core_headers.cuh"



#define MAX_EPSILON_ERROR 5e-3f



////////////////////////////////////////////////////////////////////////////////
// Constants


// Texture reference for 2D float texture
texture<float, 3, cudaReadModeElementType> tex;
const float EXTRAPVAL = 0.0f;


////////////////////////////////////////////////////////////////////////////////
//! Transform an image using texture lookups
//! @param outputData  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__global__ void transformKernel_FWD(cudaTextureObject_t thisTexObj,float *outputData,
                                uint3 dims,
                                float3 rm_1,
                                float3 rm_2,
                                float3 rm_3,
                                float3 shifts,
                                float extrapVal,
                                bool doFwdXform)
{


  // calculate normalized texture coordinates
  unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
  if (x >= dims.x) { return ; }
  unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
  if (y >= dims.y) { return ; }
  unsigned int z = blockIdx.z;
  if (z >= dims.z) { return ; }

  float u,v,w,tu,tv,tw;

  if (doFwdXform)
  {
    u = (float)x - (float)dims.x/2; 
    v = (float)y - (float)dims.y/2; 
    w = (float)z - (float)dims.z/2;
    tu = u*rm_1.x + v*rm_1.y + w*rm_1.z + shifts.x; 
    tv = u*rm_2.x + v*rm_2.y + w*rm_2.z + shifts.y;
    tw = u*rm_3.x + v*rm_3.y + w*rm_3.z + shifts.z;
  }
  else  
  {
    u = (float)x - (float)dims.x/2 - shifts.x; 
    v = (float)y - (float)dims.y/2 - shifts.y; 
    w = (float)z - (float)dims.z/2 - shifts.z;
    tu = u*rm_1.x + v*rm_1.y + w*rm_1.z; 
    tv = u*rm_2.x + v*rm_2.y + w*rm_2.z;
    tw = u*rm_3.x + v*rm_3.y + w*rm_3.z;
  }


  tu /= (float)dims.x; 
  tv /= (float)dims.y; 
  tw /= (float)dims.z;
  tu += 0.5f;
  tv += 0.5f;
  tw += 0.5f;

  if (tu < 0 | tv < 0 | tw < 0 | tu > 1 | tv > 1 | tw > 1)
  {

    outputData[ (z*dims.y + y) * dims.x + x ] = extrapVal;
  }
  else
  {
    outputData[ (z*dims.y + y) * dims.x + x ] = tex3D<float>(thisTexObj, tu, tv, tw);
  }

}


////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{

  float * d_output_img = NULL;
  float * d_input_img  = NULL;
  uint3 dims;
  float3 shifts;
  float3 rm_1;
  float3 rm_2;
  float3 rm_3;
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
    // transpose matrix
    shifts = make_float3(-ts[0],-ts[1],-ts[2]);
    rm_1   = make_float3(angles[0],angles[1],angles[2]);
    rm_2   = make_float3(angles[3],angles[4],angles[5]);
    rm_3   = make_float3(angles[6],angles[7],angles[8]);
    
  }
  else
  {

    shifts = make_float3(ts[0],ts[1],ts[2]);
    rm_1   = make_float3(angles[0],angles[3],angles[6]);
    rm_2   = make_float3(angles[1],angles[4],angles[7]);
    rm_3   = make_float3(angles[2],angles[5],angles[8]);
  }

  d_input_img = (float *)(mxGPUGetDataReadOnly(inputArray));

  mwSize const   input_dims = mxGPUGetNumberOfDimensions(inputArray);
  mwSize const * input_size = mxGPUGetDimensions(inputArray);
  mwSize const   numel_input = mxGPUGetNumberOfElements(inputArray);
  mxComplexity   input_data_type = mxGPUGetComplexity(inputArray);

  dims = make_uint3(input_size[0],input_size[1],input_size[2]);
  

  size = dims.z * dims.x * dims.y * sizeof(float);
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
  // TODO where does the 32 come from, number of bits?
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
//        cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

  checkCudaErrors(cudaMalloc3DArray(&cuArray,
                                    &channelDesc,
                                    make_cudaExtent(dims.x,dims.y,dims.z)));

  cudaMemcpy3DParms p = {0};
  p.extent = make_cudaExtent(dims.x,dims.y,dims.z);
  p.srcPtr = make_cudaPitchedPtr(d_input_img,dims.x*sizeof(float),dims.x,dims.y);
  p.dstArray = cuArray;
  p.kind = cudaMemcpyDeviceToDevice;

  cudaMemcpy3D(&p);

  struct cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(cudaResourceDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = cuArray;

  struct cudaTextureDesc texDesc;
  memset(&texDesc,0,sizeof(cudaTextureDesc));
  texDesc.filterMode = cudaFilterModeLinear;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = true;
  texDesc.addressMode[0] = cudaAddressModeWrap;
  texDesc.addressMode[1] = cudaAddressModeWrap;
  texDesc.addressMode[2] = cudaAddressModeWrap;

  cudaTextureObject_t thisTexObj;
  cudaCreateTextureObject(&thisTexObj,&resDesc,&texDesc,NULL);

//  // Set texture parameters
//  // cudaAddressModeWrap cudaAddressModeClamp cudaAddressModeMirror cudaAddressModeBorder
//  tex.addressMode[0] = cudaAddressModeWrap;
//  tex.addressMode[1] = cudaAddressModeWrap;
//  tex.addressMode[2] = cudaAddressModeWrap;

//  tex.filterMode = cudaFilterModeLinear;

//  tex.normalized = true;    // access with normalized texture coordinates



  
  // Bind the array to the texture

//  checkCudaErrors(cudaBindTextureToArray(tex3d, cuArray, channelDesc));


  // How to choose these?
//  dim3 dimBlock(8, 8, 1);
  int dimDist = 16;
  dim3 dimBlock = dim3(dimDist,dimDist,1);
  dim3 dimGrid = dim3((dims.x+dimDist-1) / dimDist, (dims.y+dimDist-1)/dimDist, dims.z);

//  dim3 dimGrid(dims.x / dimBlock.x, dims.y / dimBlock.y, 1);

  


  transformKernel_FWD<<<dimGrid, dimBlock>>>(thisTexObj, d_output_img, dims,
                                              rm_1,rm_2,rm_3,shifts,EXTRAPVAL,*doFwdXform);
  




  checkCudaErrors(cudaDeviceSynchronize());


  // Check if kernel execution generated an error
  getLastCudaError("Kernel execution failed");


  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);

  mxGPUDestroyGPUArray(inputArray);
  mxGPUDestroyGPUArray(outputArray);

  cudaFreeArray(cuArray);
  cudaDestroyTextureObject(tex3d);
//  checkCudaErrors(cudaFree(d_input_img));
//  checkCudaErrors(cudaFree(d_output_img));


}
