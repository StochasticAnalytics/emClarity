//% Matlab side code


#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cufft.h>
#include <typeinfo>


// System includes
#include <stdio.h>
#include <assert.h>

// CUDA runtime
#include <cuda.h>
#include <cuda_runtime.h>



// bool halfGrid
// bool square
// int nX,nY,nZ
// float pixelSize Ang
// float CS Ang
// float WL Ang
// float df1 Ang
// float df2 Ang
// float Ang Astig
// float AmpContrast

static __global__ void RealScale(cufftReal*, float numel_output, int hX, int nY,
                                 int oY, int oX, float WL, float CS, 
                                 float df1, float df2,float defA,float ampCont, 
                                 float fourierVoxelSize, bool doHalfGrid, bool doSqCTF) ;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{


  // Initialize the MathWorks GPU API.


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


  int oX;
  int oY;
  int hX;
  float ampCont;
  float fourierVoxelSize = 1.0f/(*pixelSize * (float)*nX);


  ampCont = atanf(*AmpContrast / sqrtf(1.0 - powf(*AmpContrast, 2)));


  cufftReal *pOut;
  // Calculate the origin for each dimension
  if (*doHalfGrid == 1)
  {
    oX = 0;
    hX = *nX/2;
  }
  else
  {
    oX = *nX/2; // int division
    hX = *nX;
  }

  oY = *nY/2;

    
  mxGPUArray *outputArray;  

  // Get the size to use in making the plan
  mwSize const   input_dims = 2;
  mxComplexity   output_type = mxREAL;
  mxClassID      output_class = mxSINGLE_CLASS;
 
  mwSize  output_size[input_dims];
  output_size[0] = (mwSize)  hX + 1;
  output_size[1] = (mwSize) *nY;


  outputArray = mxGPUCreateGPUArray(input_dims,
                                    output_size,
                                    output_class,
                                    output_type,
                                    MX_GPU_DO_NOT_INITIALIZE);

  mwSize const numel_output = mxGPUGetNumberOfElements(outputArray);

  pOut = (cufftReal *)(mxGPUGetData(outputArray));

  ////////////////////
  RealScale<<<512,256>>>(pOut, numel_output, hX + 1, *nY, oY, oX, *waveLength, *CS, 
                        *defocus1, *defocus2, *defocusAst ,ampCont, fourierVoxelSize,
                        *doHalfGrid, *doSqCTF  );


  plhs[0] = mxGPUCreateMxArrayOnGPU(outputArray);


  mxGPUDestroyGPUArray(outputArray);

}

// Complex scale
static __global__ void RealScale(cufftReal* a, float n_elements, int nX, int nY,
                                 int oY, int oX, float WL, float CS, 
                                 float defocus1, float defocus2,
                                 float defA, float ampCont, float fourierVoxelSize,
                                 bool doHalfGrid,bool doSqCTF) 
{



  // Other dimensions?   
  const int numThreads = blockDim.x * gridDim.x;
  const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  const float t1 = 3.14159265358979 * 0.5 * CS * powf(WL,3);
  const float t2 = 3.14159265358979 * 0.5 * WL;
  float radius_sq;
  float iX;
  float iY;
  int nXnY = nX*nY;
  float phi = 0;
  float phi0 = 0;

  // For now isotropic
  float df1 = 0.5*(defocus1+defocus2);
  float df2 = 0.5*(defocus1-defocus2);

//  const float scaleBy = 1/n_elements;

  if (doSqCTF == 0)
    for (int i = threadID; i < n_elements; i += numThreads) 
    {
        iY = i % (nXnY) / nX;
        iX = i - (nX)*iY;
    
        // For negative frequencies
        if (iY > oY) { iY = iY - nY ; }

        radius_sq = powf((float)iX*fourierVoxelSize,2) + powf((float)iY*fourierVoxelSize,2);
        
        a[i] = -sinf(t1*powf(radius_sq,2) + t2*(df1 + df2 * cos(2 * (phi-phi0)))* radius_sq + ampCont);
    }
  else
  {
    for (int i = threadID; i < n_elements; i += numThreads) 
    {
        iY = i % (nXnY) / nX;
        iX = i - (nX)*iY;
    
        // For negative frequencies
        if (iY > oY) { iY = iY - nY ; }
        if (doHalfGrid == 0 && iX > oX) {iX = iX - nX ;} 

        radius_sq = powf((float)iX*fourierVoxelSize,2) + powf((float)iY*fourierVoxelSize,2);
        
        a[i] = powf(sinf(t1*powf(radius_sq,2) + t2*(df1 + df2 * cos(2 * (phi-phi0)))* radius_sq + ampCont),2);
    }
  }


}

//static __global__  void AddThis(float defocus1, float defocus2, float out)
//{
//  out[0] = (defocus1[0] + defocus2[0]);
//}


