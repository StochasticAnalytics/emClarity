//% Matlab side code

#include "include/core_headers.cuh"


static __global__ void RealScale(cufftReal*, float ) ;

//static void cleanUpMemory(void);
//{
//      mexPrintf("Destroying the plans\n");
//      cufftDestroy(*plan);
//      cufftDestroy(*planInv);
//      mxGPUDestroyGPUArray(inputArray); 
//}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{



///* Check for proper number of arguments. */
//  if ( nrhs != 2) {
//      mexErrMsgIdAndTxt("MATLAB:matrixDivide:rhs",
//          "This function requires 2 input matrices.");
//  }

  // Pointers to pass to cufft
  cufftReal *pReal;
  cufftComplex *pComplex;
  mxGPUArray *outputArray;  
  int *invTrim;
  // Handle for the plan, and check to see if it exists
  cufftHandle *plan;
  cufftHandle *planInv;
  bool make_plan = true;
  bool fwd_xform = true; 
  bool do_scale = false; // Currently just as fast to do it in matlab.

/* Check for proper number of arguments. */
  if ( ! mxGPUIsValidGPUData(prhs[0]) ) 
  {

      mexErrMsgIdAndTxt("MATLAB:mexFFT:rhs",
          "This inputArray is not valid gpudata.");
  }

  // Input array, could also just pass the dimensions
  mxGPUArray const * inputArray  = mxGPUCreateFromMxArray(prhs[0]);

//  mxGPUArray const * mex_EO = mxGPUCreateFromMxArray(prhs[1]);
  invTrim =  (int *) mxGetData(prhs[1]);

  // Get the size to use in making the plan
  mwSize const   input_dims = mxGPUGetNumberOfDimensions(inputArray);
  mwSize const * input_size = mxGPUGetDimensions(inputArray);
  mwSize const   numel_input = mxGPUGetNumberOfElements(inputArray);
  mxComplexity   input_type = mxGPUGetComplexity(inputArray);

  // Assuming we are only doing either R2C or C2R
  mxComplexity output_data_type;
  if ( input_type == mxREAL ) {
//    mexPrintf("It is real bitches\n");
    output_data_type =  mxCOMPLEX ;
    pReal = (cufftReal *)(mxGPUGetDataReadOnly(inputArray));
  }
  else if (input_type == mxCOMPLEX) {
    fwd_xform = false;

    output_data_type = mxREAL ;
    pComplex = (cufftComplex *)(mxGPUGetDataReadOnly(inputArray));
  }
  else
  {
    mexErrMsgIdAndTxt("MATLAB:mexFFT:rhs",
          "This inputArray is not real or complex.");
  }

  if (nrhs > 2)
  {
    make_plan = false;


    plan    = (cufftHandle*) mxGetData(prhs[2]); 
    planInv = (cufftHandle*) mxGetData(prhs[3]); 

    if (numel_input == 1)
    {
      size_t dummy;
//      mexPrintf("Destroying the plans\n");
      if (cufftGetSize(*plan, &dummy) != CUFFT_INVALID_PLAN)
        cufftDestroy(*plan);
      // else
      //   mexPrintf("The fwd plan is not valid in destructor\n");
      if (cufftGetSize(*planInv, &dummy) != CUFFT_INVALID_PLAN)
        cufftDestroy(*planInv);
      // else
      //   mexPrintf("The inv plan is not valid in destructor\n");
      if (mxGPUIsValidGPUData(prhs[0]) )
        mxGPUDestroyGPUArray(inputArray);
      // else
      //   mexPrintf("The inputArray is not valid in destructor\n");
      return;
    }
  } 
  
  // This is only used in the inverse xform. Probably a better way. is this safe?
  // This also should probably be calculated just once in the fftTransformer class and passed in.
 
  // All transforms are done out of place, so no FFTW pafft_dimsing is assumed.
  mwSize  output_size[input_dims];
  if (input_dims > 2)
  {
    if (input_size[2] > 1)
    {
      // 3d xform 
      if (fwd_xform) { output_size[0] = input_size[0]/2+1; }
      else { output_size[0] = input_size[0]*2 - *invTrim; } 
      output_size[1] = input_size[1];
      output_size[2] = input_size[2];   
    }
    else
    {
    // 2d xform
      if (fwd_xform) { output_size[0] = input_size[0]/2+1; }
      else { output_size[0] = input_size[0]*2 - *invTrim; }
      output_size[1] = input_size[1];
      output_size[2] = 1;
    }
  }
  else if (input_dims > 1)
  {
    if (input_size[1] > 1)
      {
      // also 2d
      if (fwd_xform) { output_size[0] = input_size[0]/2+1; }
      else { output_size[0] = input_size[0]*2 - *invTrim; }
      output_size[1] = input_size[1];     
      }
    else
    {
      if (fwd_xform) { output_size[0] = input_size[0]/2+1; }
      else { output_size[0] = input_size[0]*2 - *invTrim; }
      output_size[1] = 1;   
    }
  }




  // int xFormRank;
  int fft_dims[input_dims];
  // int batchSize;
  
  if (input_dims > 2) 
  { fft_dims[2] = (int) input_size[0];
    fft_dims[1] = (int) input_size[1];
    fft_dims[0] = (int) input_size[2];
  }
  else
  {
    fft_dims[1] = (int) input_size[0];
    fft_dims[0] = (int) input_size[1];
  }

  // if (input_dims > 2) 
  // {
  //   if (fft_dims[2] > 1) { xFormRank = 3; batchSize = fft_dims[2]; } else { xFormRank = 2; batchSize = fft_dims[1]; }
  // }
  // else
  // {
  //   if (fft_dims[1] > 1) { xFormRank = 2; batchSize = fft_dims[1]; } else { xFormRank = 1; batchSize = 1;}
  // }

  outputArray = mxGPUCreateGPUArray(input_dims,
                              output_size,
                              mxGPUGetClassID(inputArray),
                              output_data_type,
                              MX_GPU_DO_NOT_INITIALIZE);

  if (fwd_xform)
  {
    pComplex = (cufftComplex *)(mxGPUGetData(outputArray));
  }
  else
  {
    pReal = (cufftReal *)(mxGPUGetData(outputArray));
  }


  mwSize const numel_output = mxGPUGetNumberOfElements(outputArray);



  if (make_plan)
  {

    mwSize const ptr_dims = 1;
    mwSize ptr_size[1];
    ptr_size[0] = (mwSize) 1;
    mxClassID output_data_class = {mxUINT32_CLASS};
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

    plan    = (cufftHandle*)mxGetData(plhs[1]);
    planInv = (cufftHandle*)mxGetData(plhs[2]);

    // Make the arrays persistent
    mexMakeArrayPersistent(plhs[1]);
    mexMakeArrayPersistent(plhs[2]);


    if ( ! fwd_xform ) {
      mexErrMsgIdAndTxt("The cufft plan must be created on a forward transform.","womp");
    }

    if (input_dims == 3) {
       cufftPlan3d(plan, fft_dims[0], fft_dims[1], fft_dims[2],  CUFFT_R2C);
       cufftPlan3d(planInv, fft_dims[0], fft_dims[1],fft_dims[2], CUFFT_C2R);
    }
    else if (input_dims == 2) {
       cufftPlan2d(plan, fft_dims[0], fft_dims[1],  CUFFT_R2C);
       cufftPlan2d(planInv, fft_dims[0], fft_dims[1], CUFFT_C2R);
    }
    else {
      mexErrMsgIdAndTxt("The cufft plan must be created on a 1, 2, or 3 dimensional array.","womp");
    }


    // cufftPlanMany(plan,    xFormRank, fft_dims, 
    //               NULL, NULL, NULL, NULL, NULL, NULL,
    //               CUFFT_R2C, 1);
    // cufftPlanMany(planInv, xFormRank, fft_dims, 
    //               NULL, NULL, NULL, NULL, NULL, NULL,
    //               CUFFT_C2R, 1);

  }

  // Do the fft
  if (fwd_xform)
  {
//    mexPrintf("Doing the forward xform\n");
    cufftExecR2C(*plan, (cufftReal *)pReal, (cufftComplex *)pComplex);
  }
  else
  {
//    mexPrintf("Doing the inverse xform\n");
    cufftExecC2R(*planInv, pComplex, pReal);
    // Should probably be an option, but enforce scaling on the inverse
    // I don't know how to properly chose the block/thread configuration.
    if (do_scale)
    {
      RealScale<<<512,256>>>(pReal, numel_output);
    }
  }

  // mxGPUCreateMxArrayOnGPU is used to wrap the output data in an MxArray that can be passed back to matlab.
  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);


//  mexPrintf("Afft_dimsress of plan is %d\n", *plan);
//  mexPrintf("Afft_dimsress of planInv is %d\n", *planInv);
  mxGPUDestroyGPUArray(inputArray);
  mxGPUDestroyGPUArray(outputArray);
//  mxGPUDestroyGPUArray(mex_EO);
}

// Complex scale
static __global__ void RealScale(cufftReal *a,  float n_elements) 
{



  // Other dimensions?   
  const int numThreads = blockDim.x * gridDim.x;
  const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  const float scaleBy = 1/n_elements;

  for (int i = threadID; i < n_elements; i += numThreads) {
    a[i] *= scaleBy;
  }



}


