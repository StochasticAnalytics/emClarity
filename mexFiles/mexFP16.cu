
#include "include/core_headers.cuh"

// #define mexFP16_DEBUG_PRINT(args) mexPrintf("%s\n", args)
#define mexFP16_DEBUG_PRINT(...)

__global__ void convert_fp16_to_fp32(const uint16_t* __restrict__ input_half, float* __restrict__ output_single, const int N) {
  // Could be improved with a simple vector load.
  for (int idx = blockIdx.x * blockDim.x + threadIdx.x;  idx < N; idx += gridDim.x * blockDim.x) 
    output_single[idx] = __half2float(__ushort_as_half(input_half[idx]));
  
}

__global__ void convert_fp32_to_fp16(const float* __restrict__ input_single, uint16_t* __restrict__ output_half, const int N) {
  // Could be improved with a simple vector store.
  for (int idx = blockIdx.x * blockDim.x + threadIdx.x;  idx < N; idx += gridDim.x * blockDim.x) 
     output_half[idx] = __half_as_ushort(__float2half_rn(input_single[idx]));
  
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {

  if (nrhs != 4) {
    mexErrMsgIdAndTxt("MATLAB:mexFP16:rhs",
                      "This function requires 2 input matrices, a boolean (to_half), and an int (*n_elements).");
  }
  bool single_array_is_on_gpu = false;
  bool half_array_is_on_gpu = false;
  // First check to see if we have a gpu arra
  if (mxIsGPUArray(prhs[0])) {
    // Now let's see if it isvalid data
    if (!mxGPUIsValidGPUData(prhs[0])) {
      mexErrMsgIdAndTxt("MATLAB:mexFP16:rhs",
                        "This function requires the first input to be a valid gpuArray.");
    }
    // And if it is single precision ( written this way causes a segfault, trust that emc_halfcast has done the checking.)
    // if (mxGPUGetClassID((const mxGPUArray *)prhs[0]) != mxSINGLE_CLASS) {
    //   mexErrMsgIdAndTxt("MATLAB:mexFP16:rhs",
    //                     "This function requires the first input to be of type single mxGPUGetClassID.");
    // }
    single_array_is_on_gpu = true;  
  }
  
  // Same thing for the half array
  if (mxIsGPUArray(prhs[1])) {
    if (!mxGPUIsValidGPUData(prhs[1])) {
      mexErrMsgIdAndTxt("MATLAB:mexFP16:rhs",
                        "This function requires the second input to be a valid gpuArray.");
    }
    //( written this way causes a segfault, trust that emc_halfcast has done the checking.)
    // if (mxGPUGetClassID((const mxGPUArray *)prhs[1]) != mxUINT16_CLASS) {
    //    mexPrintf("Here %d\n", __LINE__);
    //   mexErrMsgIdAndTxt("MATLAB:mexFP16:rhs",
    //                     "This function requires the second input to be of type uint16 mxGPUGetClassID.");
    // }
    half_array_is_on_gpu = true;
  }

  if (!mxIsLogical(prhs[2])) {
    mexErrMsgIdAndTxt("MATLAB:mexFP16:rhs",
                      "This function requires the third input to be of type logical.");
  }

  if (!mxIsInt64(prhs[3])) {
    mexErrMsgIdAndTxt("MATLAB:mexFP16:rhs",
                      "This function requires the fourth input to be of type int.");
  }


  float* input_single;
  uint16_t* input_uint16;

  // mxGPUCreateFromMxArray will return a read only pointer if the underlying is a matlab gpuArray
  // To avoid a copy but still get a pointer we can cast the const away.
  // This seems dodgy as fuck
  if (single_array_is_on_gpu) {
    mxGPUArray const * inputArray  = mxGPUCreateFromMxArray(prhs[0]);
    input_single =  (float *) mxGPUGetData((mxGPUArray *)inputArray);
  }
  else
    input_single =  (float *) mxGetData(prhs[0]);

  if (half_array_is_on_gpu) {
    mxGPUArray const * inputArray  = mxGPUCreateFromMxArray(prhs[1]);
    input_uint16 = (uint16_t *) mxGPUGetData((mxGPUArray *)inputArray);
  }
  else
    input_uint16 = (uint16_t *) mxGetData(prhs[1]);

  bool* cast_to_fp16 = (bool *) mxGetData(prhs[2]);
  size_t* n_elements = (size_t *) mxGetData(prhs[3]);



  // First, we need to take care of any data conversion needed.
  float* temporary_single;
  uint16_t* temporary_uint16;

  // Output half on gpu so we need to use the conversion kernel
  if (cast_to_fp16[0]) {
    mexFP16_DEBUG_PRINT("Casting to half\n");
    if (half_array_is_on_gpu) {
      mexFP16_DEBUG_PRINT("Casting to half that is on device already\n");
      temporary_uint16 = input_uint16;
      if (single_array_is_on_gpu) {
        temporary_single = input_single;
      } 
      else {
        mexFP16_DEBUG_PRINT("Copying single to device\n");
        checkCudaErrors(cudaMallocAsync(&temporary_single, *n_elements * sizeof(float), cudaStreamPerThread));
        checkCudaErrors(cudaMemcpyAsync(temporary_single, input_single, *n_elements  * sizeof(float), cudaMemcpyHostToDevice, cudaStreamPerThread));
      }
      convert_fp32_to_fp16<<<1, 1, 0, cudaStreamPerThread>>>(temporary_single, temporary_uint16, *n_elements);

      if (!single_array_is_on_gpu)
        checkCudaErrors(cudaFreeAsync(temporary_single, cudaStreamPerThread));

      checkCudaErrors(cudaStreamSynchronize(cudaStreamPerThread));
    }
    else {
      mexFP16_DEBUG_PRINT("Casting to half that is on host already\n");
      // Even though we could do this in place, then only copy over the half precision data to the host,
      // we want to keep the input data un altered.
      if (single_array_is_on_gpu) {
        mexFP16_DEBUG_PRINT("Copying single to host\n");
        checkCudaErrors(cudaMallocHost(&temporary_single, *n_elements * sizeof(float)));
        checkCudaErrors(cudaMemcpyAsync(temporary_single, input_single, *n_elements  * sizeof(float), cudaMemcpyDeviceToHost, cudaStreamPerThread));
        checkCudaErrors(cudaStreamSynchronize(cudaStreamPerThread));
      } 
      else 
        temporary_single = input_single;

      mexFP16_DEBUG_PRINT(std::to_string(*n_elements));
      half_float::half* half_ptr = (half_float::half *)(input_uint16);

      for (int i = 0; i < *n_elements ; i++) {
        half_ptr[i] = half_float::half(temporary_single[i]);
      }

      if (single_array_is_on_gpu)
        checkCudaErrors(cudaFreeHost(temporary_single));

    }
  }
  else {
    mexFP16_DEBUG_PRINT("Casting to single\n");
    // Casting from half to single
    if (single_array_is_on_gpu) {
      temporary_single = input_single;
      if (half_array_is_on_gpu) {
        temporary_uint16 = input_uint16;
      } 
      else {
        checkCudaErrors(cudaMallocAsync(&temporary_uint16, *n_elements * sizeof(uint16_t), cudaStreamPerThread));
        checkCudaErrors(cudaMemcpyAsync(temporary_uint16, input_uint16, *n_elements  * sizeof(uint16_t), cudaMemcpyHostToDevice, cudaStreamPerThread));
      }
      convert_fp16_to_fp32<<<1024, 32, 0, cudaStreamPerThread>>>(temporary_uint16, temporary_single, *n_elements);

      if (!half_array_is_on_gpu)
        checkCudaErrors(cudaFreeAsync(temporary_uint16, cudaStreamPerThread));

     checkCudaErrors(cudaStreamSynchronize(cudaStreamPerThread));
    }
    else {
      mexFP16_DEBUG_PRINT("Casting to single that is on host already\n");
      // the output single is on the host
      half_float::half* half_ptr = (half_float::half *)(input_uint16);
      if (half_array_is_on_gpu) {
        mexFP16_DEBUG_PRINT("Copying half to host\n");
        checkCudaErrors(cudaMallocHost(&temporary_uint16, *n_elements * sizeof(uint16_t)));
        checkCudaErrors(cudaMemcpyAsync(temporary_uint16, input_uint16, *n_elements  * sizeof(uint16_t), cudaMemcpyDeviceToHost, cudaStreamPerThread));
        checkCudaErrors(cudaStreamSynchronize(cudaStreamPerThread));
        half_ptr = (half_float::half *)(temporary_uint16);
      } 

      for (int i = 0; i < *n_elements ; i++) {
        float tmp = half_float::half_cast<float>(half_ptr[i]);
        input_single[i] = half_float::half_cast<float>(half_ptr[i]);
      }

      if (half_array_is_on_gpu)
        checkCudaErrors(cudaFreeHost(temporary_uint16));
    }
  } // end if if cast to fp16 else to single


  return;

    
}