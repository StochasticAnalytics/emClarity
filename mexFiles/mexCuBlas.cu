#include "include/core_headers.cuh"



#define M 2
#define N 3
#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{

  // create or destroy a cublas handle

  RotationMatrix rotmat;
  RotationMatrix updateMat;
  rotmat.SetToEulerRotation(1,20,30);
  updateMat.SetToEulerRotation(30,40,-20);
  mexPrintf("\nForward\n%2.4f %2.4f %2.4f\n%2.4f %2.4f %2.4f\n%2.4f %2.4f %2.4f\n\n",rotmat.m[0][0],rotmat.m[0][1],rotmat.m[0][2],rotmat.m[1][0],rotmat.m[1][1],rotmat.m[1][2],rotmat.m[2][0],rotmat.m[2][1],rotmat.m[2][2]);

  rotmat *= updateMat;


  mexPrintf("\nUpdated\n%2.4f %2.4f %2.4f\n%2.4f %2.4f %2.4f\n%2.4f %2.4f %2.4f\n\n",rotmat.m[0][0],rotmat.m[0][1],rotmat.m[0][2],rotmat.m[1][0],rotmat.m[1][1],rotmat.m[1][2],rotmat.m[2][0],rotmat.m[2][1],rotmat.m[2][2]);

  bool createHandle = false;
  cublasHandle_t *handle; 
  if ( nrhs == 0 )
  { 
    mexPrintf("Creating a cublas handle\n"); 
    createHandle = true;
  }

  // Initialize the MathWorks GPU API.
  mxInitGPU();



    cudaError_t cudaStat;    
    cublasStatus_t stat;


    int i, j;
    float* devPtrA;
    float* a = 0;
    float* sumVal;
    a = (float *)malloc (M * N * sizeof (*a));
    if (!a) {
        printf ("host memory allocation failed");
        return;
    }
    for (j = 1; j <= N; j++) {
        for (i = 1; i <= M; i++) {
            a[IDX2F(i,j,M)] = (float)((i-1) * M + j);
            mexPrintf("%d %d %f\n",i,j,a[IDX2F(i,j,M)]);
        }
    }
    cudaStat = cudaMalloc ((void**)&devPtrA, M*N*sizeof(*a));
    if (cudaStat != cudaSuccess) {
        printf ("device memory allocation failed");
        return;
    }

    // Create the cublas context if needed
    if (createHandle) 
    { 

      mwSize const ptr_dims = 1;
      mwSize ptr_size[1];
      ptr_size[0] = (mwSize) 1;
      mxClassID output_data_class = {mxUINT32_CLASS};
      mxComplexity output_data_complexity = {mxREAL};

      // Forward and inverse transforms
      plhs[0] =  mxCreateNumericArray(ptr_dims,
                                      ptr_size,
                                      output_data_class,
                                      output_data_complexity);


      handle    = (cublasHandle_t *)mxGetData(plhs[0]);

      // Make the arrays persistent
      mexMakeArrayPersistent(plhs[0]);
      stat = cublasCreate(handle); 
mexPrintf("%p\n" , handle);
    }
    else
    {
      mexPrintf("Reusing the handle\n");
      handle    = (cublasHandle_t *) mxGetData(prhs[0]); 
mexPrintf("%p\n" , handle);
    }


    stat = cublasSetMatrix (M, N, sizeof(*a), a, M, devPtrA, M);
    
    float alpha = 3.0f;
    cublasSscal( *handle, M*N, &alpha, devPtrA, 1);
//    for (j = 1; j <= N; j++) {
//        for (i = 1; i <= M; i++) {
//            printf ("%7.0f", a[IDX2F(i,j,M)]);
//        }
//        printf ("\n");
//    }

    free(a);
    return;



}
