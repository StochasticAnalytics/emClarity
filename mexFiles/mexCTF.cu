//% Matlab side code

#include "include/core_headers.cuh"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{


  // Initialize the MathWorks GPU API.
  mxInitGPU();

  // TODO if these aren't typed properly in the mex call, I can't cast them appropriatley here
  bool * doHalfGrid =  (bool *) mxGetData(prhs[0]);  
  bool * doSqCTF    =  (bool *) mxGetData(prhs[1]); 
  int  * nX         =  (int  *) mxGetData(prhs[2]);
  int  * nY         =  (int  *) mxGetData(prhs[3]);
  float* pixelSize  =  (float*) mxGetData(prhs[4]);
  float* waveLength =  (float*) mxGetData(prhs[5]);
  float* CS         =  (float*) mxGetData(prhs[6]);
  float* defocus1   =  (float*) mxGetData(prhs[7]);
  float* defocus2   =  (float*) mxGetData(prhs[8]);
  float* defocusAst =  (float*) mxGetData(prhs[9]);
  float* AmpContrast=  (float*) mxGetData(prhs[10]);
  bool calc_centered = false;
  float* occupancy = NULL;
  float* exposure = NULL;
  bool use_weights = false;

  if ( nrhs == 12) 
  { 
    mexPrintf("Doing centered ctf calc\n");
    bool * tmp_bool = (bool *) mxGetData(prhs[11]);
    calc_centered = *tmp_bool;
  }

  if ( nrhs == 14) 
  { 
    bool * tmp_bool = (bool *) mxGetData(prhs[11]);
    calc_centered = *tmp_bool;
    use_weights = true;
    occupancy = (float *) mxGetData(prhs[12]);
    exposure = (float *) mxGetData(prhs[13]);
  }


  ctfParams b_ctf(*doHalfGrid,*doSqCTF,*pixelSize,*waveLength,*CS,*AmpContrast,
                  *defocus1,  *defocus2, *defocusAst);

//  mexPrintf("%f %f %f %f\n",*defocus1,*defocus2,b_ctf.defocus1,b_ctf.defocus2);
  uint2 dims;
  uint2 o_dims;

  

  float2 fourierVoxelSize;
  fourierVoxelSize = make_float2( 1.0f/(*pixelSize * (float)*nX), 1.0f/(*pixelSize * (float)*nY));
 


  cufftReal *pOut;
  // Calculate the origin for each dimension
  // Get the size to use in making the plan
  mwSize const   input_dims = 2;
  mxComplexity   output_type = mxREAL;
  mxClassID      output_class = mxSINGLE_CLASS;
  dims = make_uint2(*nX,*nY);
  mwSize  output_size[input_dims];

  if (*doHalfGrid )
  {
    o_dims = make_uint2(0, *nY/2);
    dims.x = dims.x/2 + 1;
  }
  else
  {
    o_dims = make_uint2(*nX/2, *nY/2);

  }

  output_size[0] = (mwSize)  dims.x;
  output_size[1] = (mwSize)  dims.y;
    
  mxGPUArray *outputArray;  
  // TODO it would be nice not to init all the zeros, but then the fourier padding would need to be dealt with.
  outputArray = mxGPUCreateGPUArray(input_dims,
                                    output_size,
                                    output_class,
                                    output_type,
                                    MX_GPU_INITIALIZE_VALUES);


  pOut = (float *)(mxGPUGetData(outputArray));

  ////////////////////
  dim3 dimBlock(32, 32, 1);
  dim3 dimGrid((dims.x + dimBlock.x -1) / dimBlock.x, (dims.y + dimBlock.y -1) / dimBlock.y, 1);

//  mexPrintf("%d %d %d %d\n",dims.x, dims.y, o_dims.x, o_dims.y);
//  mexPrintf("%f %f\n",fourierVoxelSize.x, fourierVoxelSize.y);
//  const  float t1 =  PI * 0.5f * b_ctf.CS* powf(b_ctf.waveLength,3);
//  const  float t2  = PI * b_ctf.waveLength;
//  mexPrintf("%f %f\n",t1,t2);

//  if (use_weights)
//  {
//    ctf<<< dimGrid, dimBlock >>>(pOut, dims, o_dims, b_ctf, fourierVoxelSize, calc_centered,*occupancy,*exposure);
//  }
//  else
//  {
    ctf<<< dimGrid, dimBlock >>>(pOut, dims, o_dims, b_ctf, fourierVoxelSize, calc_centered);
//  }



  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);


  mxGPUDestroyGPUArray(outputArray);

}


