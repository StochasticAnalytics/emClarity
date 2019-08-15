// Mex includes
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "matrix.h"

// Cuda includes
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "cublas_v2.h"

// Utilities and timing functions
#include "helper_functions.h"    // includes cuda.h and cuda_runtime_api.h
#include "helper_cuda.h"         // helper functions for CUDA error check

#include <typeinfo>

// System includes
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "../utils/rotation_matrix.h"



const float PI = 3.14159265358979323;
const float PI_sq = PI * PI;
const float PI_half = 0.5f*PI;

inline float deg_2_rad(float input_angle) { return input_angle * PI / 180.0f ;}
inline float rad_2_deg(float input_angle) { return input_angle * 180.0f / PI ;}


// A drop in for debugging call with c++ macro __LINE__ as the argument
inline void mexit(int line_number){ mexPrintf("\n\tAt line %d\n\n", line_number); mexErrMsgIdAndTxt("MATLAB:mexSF3D", "Dropping out per mexit") ; }


struct ctfParams {

  bool  doHalfGrid;
  bool  doSqCTF;
  float pixelSize; // Angstrom

  float waveLength; // Angstrom
  float CS; // millimeter
  float amplitudeContrast;
  float defocus1; // Angstrom
  float defocus2; // Angstrom
  float astigmatism_angle; // from x-axis

	  __host__ __device__ ctfParams() : doHalfGrid(true), doSqCTF(false),
                                      pixelSize(0.0f), waveLength(0.0f), 
                                      CS(0.0f), amplitudeContrast(0.0f), 
                                      defocus1(0.0f), defocus2(0.0f), 
                                      astigmatism_angle(0.0f) {}
	  __host__ __device__ ctfParams(bool doHalfGrid, bool doSqCTF,
                                  float pixelSize, float waveLength,
                                  float CS, float amplitudeContrast, 
                                  float df1, float df2, 
                                  float astigmatism_angle) : 
                                  doHalfGrid(doHalfGrid), doSqCTF(doSqCTF),
                                  pixelSize(pixelSize), waveLength(waveLength),
                                  CS(CS * 1e7), 
                                  amplitudeContrast(atanf(amplitudeContrast / 
                                                          sqrtf(1.0 - powf(amplitudeContrast, 2)))), // Convert ampContrast to phase shift TODO is this safe? 
                                  defocus1(0.5f*(df1+df2)), defocus2(0.5f*(df1-df2)), // Convert these to terms used in calc
                                  astigmatism_angle(astigmatism_angle*PI/180.0f - (PI * (float)lrintf(astigmatism_angle/180.0f))) {} // enforce -90 to 90 convention





};


// Kernel defs

__global__ void ctf(cufftReal* a, uint2 dims, uint2 o_dims, ctfParams b_ctf, float2 fourierVoxelSize, bool calc_centered);
