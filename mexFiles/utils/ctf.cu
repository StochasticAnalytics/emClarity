// CUDA runtime
//#include "../include/core_headers.cuh"
#include "../include/bh_helper.cuh"

// createb_ctf
__global__ void ctf(cufftReal* a, uint2 dims, uint2 o_dims, ctfParams b_ctf, float2 fourierVoxelSize, bool calc_centered) 
{


  // Maybe put this in the ctf struct
  const  float t1 =  PI * 0.5f * b_ctf.CS* powf(b_ctf.waveLength,3);
  const  float t2  = PI * b_ctf.waveLength;

  // Other dimensions?   
  int x = blockIdx.x*blockDim.x + threadIdx.x;
  int y = blockIdx.y*blockDim.y + threadIdx.y;

  float tmp_coord;
  float radius_sq;
  long output_IDX;
  float phi;

  // TODO add centered calc
  output_IDX = y*dims.x + x;

  if (calc_centered)
  {
    y -= o_dims.y;
    if ( ! b_ctf.doHalfGrid ) { x -= o_dims.x ;}
  }
  else
  {
    if (y > o_dims.y) { y = y - dims.y ; }
    if ( ! b_ctf.doHalfGrid  && x > o_dims.x) {x = x - dims.x ;}
  }
 


  // TODO this seems safe, buy could there be probs for x or y == 0??
  tmp_coord = (float)y*fourierVoxelSize.y;
  radius_sq = (float)x*fourierVoxelSize.x;
  phi = atan2f(tmp_coord,radius_sq); 

  radius_sq = radius_sq*radius_sq + tmp_coord*tmp_coord;
  
  a[output_IDX] = sinf(t1*powf(radius_sq,2) + t2*radius_sq*(b_ctf.defocus1 + b_ctf.defocus2 * cosf(2.0f * (phi-b_ctf.astigmatism_angle))) + b_ctf.amplitudeContrast);

  if (b_ctf.doSqCTF )
  {
    // Is this any better (or worse) than pow?
    a[output_IDX] *= a[output_IDX];
  }
    
}

