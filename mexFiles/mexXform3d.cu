#include "include/core_headers.cuh"



#define MAX_EPSILON_ERROR 5e-3f



////////////////////////////////////////////////////////////////////////////////
// Constants


// Texture reference for 2D float texture
//texture<float, 3, cudaReadModeElementType> tex;
const float EXTRAPVAL = 0.0f;


////////////////////////////////////////////////////////////////////////////////
//! Transform an image using texture lookups
//! @param outputData  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__global__ void transformKernel_FWD(cudaTextureObject_t thisTexObj,
                                    float* outputData,
                                    int3 dims,
                                    float3 rm_1,
                                    float3 rm_2,
                                    float3 rm_3,
                                    float3 shifts,
                                    float extrapVal,
                                    bool doFwdXform,
                                    float3 size_shift)
{


  // calculate normalized texture coordinates
  int x = blockIdx.x*blockDim.x + threadIdx.x;
  if (x >= dims.x) { return ; }
  int y = blockIdx.y*blockDim.y + threadIdx.y;
  if (y >= dims.y) { return ; }
  int z = blockIdx.z*blockDim.z + threadIdx.z;
  if (z >= dims.z) { return ; }

  float u,v,w,tu,tv,tw;

  if (doFwdXform)
  {

    // First, put the origin at the center, rotate, then shift
    u = (float)x - (float)(dims.x/2) + shifts.x; 
    v = (float)y - (float)(dims.y/2) + shifts.y; 
    w = (float)z - (float)(dims.z/2) + shifts.z;
    tu = u*rm_1.x + v*rm_1.y + w*rm_1.z + size_shift.x; 
    tv = u*rm_2.x + v*rm_2.y + w*rm_2.z + size_shift.y;
    tw = u*rm_3.x + v*rm_3.y + w*rm_3.z + size_shift.z;
  }
  else  
  {
    // First, put the origin at the center, shift, then rotate
    u = (float)x - (float)(dims.x/2); 
    v = (float)y - (float)(dims.y/2); 
    w = (float)z - (float)(dims.z/2);
    tu = u*rm_1.x + v*rm_1.y + w*rm_1.z + shifts.x + size_shift.x; 
    tv = u*rm_2.x + v*rm_2.y + w*rm_2.z + shifts.y + size_shift.y;
    tw = u*rm_3.x + v*rm_3.y + w*rm_3.z + shifts.z + size_shift.z;
  }

  // Convert to normalized coordinates 
  tu /= (float)dims.x; 
  tv /= (float)dims.y; 
  tw /= (float)dims.z;
  tu += 0.5f;
  tv += 0.5f;
  tw += 0.5f;

//  if (tu < 0 | tv < 0 | tw < 0 | tu >= 1 - 1/(float)dims.x | tv >= 1 - 1/(float)dims.y | tw >= 1 - 1/(float)dims.z)
//  {

//    outputData[ (z*dims.y + y) * dims.x + x ] = extrapVal;
//  }
//  else
//  {
//    outputData[ (z*dims.y + y) * dims.x + x ] = tex3D<float>(thisTexObj, tu, tv, tw);
//  }
    // cudaAddressModeBorder returns 0.0 for out of bounds reads.
    outputData[ (z*dims.y + y) * dims.x + x ] = tex3D<float>(thisTexObj, tu, tv, tw);
}


////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{

  float * d_output_img = NULL;
  float * d_input_img  = NULL;
  int3 dims;
  float3 shifts;
  float3 rm_1;
  float3 rm_2;
  float3 rm_3;
  bool *doFwdXform; // rotate then shift



  // We want to be able to re-use the texture object, so only set it up once.
  bool make_texture_obj = true;
  cudaTextureObject_t* tex;
  cudaArray_t *cuArray;
  

/* Check for proper number of arguments. TODO add checks on narg and types*/

  if ( nrhs == 2 )
  {


    tex    = (cudaTextureObject_t *) mxGetData(prhs[0]);   
    checkCudaErrors(cudaDestroyTextureObject(*tex));

    cuArray = (cudaArray_t *) mxGetData(prhs[1]);
    checkCudaErrors(cudaFreeArray(*cuArray));

//   mxFree(tex);
    return;
  }

  if ( ! mxGPUIsValidGPUData(prhs[1])) 
  {
        mexErrMsgIdAndTxt("MATLAB:mexFFT:rhs", "Expected a 3d gpuarray.");  
  }



  if ( nrhs > 5 )
  {
      tex    = (cudaTextureObject_t *) mxGetData(prhs[5]);   
      make_texture_obj = false; 

  }

  // This should be the input array or the the original arrays dims
  mxGPUArray const * inputArray  = mxGPUCreateFromMxArray(prhs[1]);
  mwSize const   input_dims = mxGPUGetNumberOfDimensions(inputArray);
  mxComplexity   input_data_type = mxGPUGetComplexity(inputArray);
  
  size_t* input_size = (size_t *) mxGetData(prhs[0]);  


  mxGPUArray * outputArray;

  // The first arg is either a 3d or a pointer to the tex obj w/ previous 3d passed.
  float *angles = (float *) mxGetData(prhs[2]);



  float *ts = (float *) mxGetData(prhs[3]);
  doFwdXform = (bool *) mxGetData(prhs[4]);

  if (*doFwdXform)
  {
    // transpose matrix
    shifts = make_float3(-ts[0],-ts[1],-ts[2]);
  }
  else
  {
    shifts = make_float3(ts[0],ts[1],ts[2]);
  }

    rm_1   = make_float3(angles[0],angles[3],angles[6]);
    rm_2   = make_float3(angles[1],angles[4],angles[7]);
    rm_3   = make_float3(angles[2],angles[5],angles[8]);

//  mexPrintf("ts %f %f %f\n", shifts.x, shifts.y, shifts.z);
//mexPrintf("rm1 %f %f %f\n", rm_1.x, rm_1.y, rm_1.z);
//mexPrintf("rm2 %f %f %f\n", rm_2.x, rm_2.y, rm_2.z);
//mexPrintf("rm3 %f %f %f\n", rm_3.x, rm_3.y, rm_3.z);

  dims = make_int3(input_size[0],input_size[1],input_size[2]);
  



  // Create MX array and init with zeros
  outputArray = mxGPUCreateGPUArray(input_dims,
                                    input_size,
                                    mxSINGLE_CLASS,
                                    input_data_type,
                                    MX_GPU_DO_NOT_INITIALIZE);

  d_output_img = (float *)(mxGPUGetData(outputArray));


  if (make_texture_obj)
  {

    d_input_img = (float *)(mxGPUGetDataReadOnly(inputArray));

    mwSize const ptr_dims = 1;
    mwSize ptr_size[1];
    ptr_size[0] = (mwSize) 1;
    mxClassID output_data_class = {mxUINT64_CLASS};
    mxComplexity output_data_complexity = {mxREAL};

    // Forward and inverse transforms
    plhs[1] =  mxCreateNumericArray(ptr_dims,
                                    ptr_size,
                                    output_data_class,
                                    output_data_complexity);

    plhs[2] =  mxCreateNumericArray(ptr_dims,
                                    ptr_size,
                                    output_data_class,
                                    output_data_complexity);

    tex = (cudaTextureObject_t *)mxGetData(plhs[1]);

    cuArray = (cudaArray_t *)mxGetData(plhs[2]);


    *tex = 0;

    // Make the arrays persistent
//    mexMakeArrayPersistent(plhs[1]);
//    mexMakeArrayPersistent(plhs[2]);
//  

//    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

    checkCudaErrors(cudaMalloc3DArray(cuArray,
                                      &channelDesc,
                                      make_cudaExtent(dims.x,dims.y,dims.z)));

    cudaMemcpy3DParms p = {0};
    p.extent = make_cudaExtent(dims.x,dims.y,dims.z);
    p.srcPtr = make_cudaPitchedPtr(d_input_img, dims.x*sizeof(float),dims.x,dims.y);
    p.dstArray = *cuArray;
    p.kind = cudaMemcpyDeviceToDevice;

    cudaMemcpy3D(&p);

//    mexPrintf(" cuArray in allocate %p %p %p\n", cuArray, &cuArray, *cuArray);

    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = *cuArray;
   

//    size_t n_bytes = dims.x*dims.y*dims.z*sizeof(float);
//    float* buffer;
//    checkCudaErrors(cudaMalloc(&buffer, n_bytes );
//    struct cudaResourceDesc resDesc;

//    memset(&resDesc, 0, sizeof(cudaResourceDesc));
//    resDesc.resType = cudaResourceTypeLinear;
//    resDesc.res.linear.devPtr = buffer;
//    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
//    resDesc.res.linear.desc.x = 32; // bits per channel
//    resDesc.res.linear.sizeInBytes = n_bytes;

    struct cudaTextureDesc texDesc;
    memset(&texDesc,0,sizeof(texDesc));

    texDesc.filterMode = cudaFilterModeLinear;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = 1;
    texDesc.addressMode[0] = cudaAddressModeBorder; //cudaAddressModeClamp;
    texDesc.addressMode[1] = cudaAddressModeBorder; //cudaAddressModeClamp;
    texDesc.addressMode[2] = cudaAddressModeBorder; //cudaAddressModeClamp;

    


    // TODO does the cuArray need to be persistent?


    cudaCreateTextureObject(tex,&resDesc,&texDesc,NULL);



  }



  dim3 dimBlock = dim3(32,32,1);
  dim3 dimGrid = dim3((dims.x+32-1) / 32, (dims.y+32-1)/32, dims.z);


//  mexPrintf("Dims: %d %d %d\nBlock: %d %d %d\nGrid: %d %d %d\n",
//            dims.x, dims.y, dims.z, dimBlock.x, dimBlock.y, dimBlock.z, dimGrid.x, dimGrid.y, dimGrid.z);

  float3 size_shift = make_float3(0.0f, 0.0f, 0.0f);
  if (IsEven(dims.x)) size_shift.x = 0.5f;
  if (IsEven(dims.y)) size_shift.y = 0.5f;
  if (IsEven(dims.z)) size_shift.z = 0.5f;

  checkCudaErrors(cudaGetLastError());
  transformKernel_FWD<< <dimGrid, dimBlock, 0, cudaStreamDefault>> >(*tex, d_output_img, dims,
                                              rm_1,rm_2,rm_3,shifts,EXTRAPVAL,*doFwdXform, size_shift);
  cudaError_t cuda_error = cudaStreamSynchronize(cudaStreamDefault); 
  if (cuda_error != cudaSuccess) 
  {
    mexPrintf("Sync Check error = %s at line %d\n", _cudaGetErrorEnum(cuda_error), __LINE__);
  }



  // Check if kernel execution generated an error
  getLastCudaError("Kernel execution failed");


  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);

  mxGPUDestroyGPUArray(inputArray);
  mxGPUDestroyGPUArray(outputArray);

//  cudaFreeArray(cuArray);
//  cudaDestroyTextureObject(tex3d);

//  checkCudaErrors(cudaFree(d_output_img));


}
