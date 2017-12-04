#ifndef _MACROSH_
#define _MACROSH_

//Need to be sure that Nx,Ny,Nz are defined within scope to use these
//General indexing...
#define GET3DINDEX_XYZ k*Nx*Ny + j*Nx + i
#define GET3DINDEX_YZX i*Nz*Ny + k*Ny + j
#define GET3DINDEX_ZXY j*Nx*Nz + i*Nz + k

//For pulling data for Neumann BC's
#define GET3DINDEX_XYZ_Xp1 k*Nx*Ny + j*Nx + i+1
#define GET3DINDEX_XYZ_Xp2 k*Nx*Ny + j*Nx + i+2
#define GET3DINDEX_XYZ_Xp3 k*Nx*Ny + j*Nx + i+3
#define GET3DINDEX_XYZ_Xp4 k*Nx*Ny + j*Nx + i+4
#define GET3DINDEX_XYZ_Xp5 k*Nx*Ny + j*Nx + i+5
#define GET3DINDEX_XYZ_Xp6 k*Nx*Ny + j*Nx + i+6
#define GET3DINDEX_XYZ_Xm1 k*Nx*Ny + j*Nx + i-1
#define GET3DINDEX_XYZ_Xm2 k*Nx*Ny + j*Nx + i-2
#define GET3DINDEX_XYZ_Xm3 k*Nx*Ny + j*Nx + i-3
#define GET3DINDEX_XYZ_Xm4 k*Nx*Ny + j*Nx + i-4
#define GET3DINDEX_XYZ_Xm5 k*Nx*Ny + j*Nx + i-5
#define GET3DINDEX_XYZ_Xm6 k*Nx*Ny + j*Nx + i-6

#define GET3DINDEX_XYZ_Yp1 k*Nx*Ny + (j+1)*Nx + i
#define GET3DINDEX_XYZ_Yp2 k*Nx*Ny + (j+2)*Nx + i
#define GET3DINDEX_XYZ_Yp3 k*Nx*Ny + (j+3)*Nx + i
#define GET3DINDEX_XYZ_Yp4 k*Nx*Ny + (j+4)*Nx + i
#define GET3DINDEX_XYZ_Yp5 k*Nx*Ny + (j+5)*Nx + i
#define GET3DINDEX_XYZ_Yp6 k*Nx*Ny + (j+6)*Nx + i
#define GET3DINDEX_XYZ_Ym1 k*Nx*Ny + (j-1)*Nx + i
#define GET3DINDEX_XYZ_Ym2 k*Nx*Ny + (j-2)*Nx + i
#define GET3DINDEX_XYZ_Ym3 k*Nx*Ny + (j-3)*Nx + i
#define GET3DINDEX_XYZ_Ym4 k*Nx*Ny + (j-4)*Nx + i
#define GET3DINDEX_XYZ_Ym5 k*Nx*Ny + (j-5)*Nx + i
#define GET3DINDEX_XYZ_Ym6 k*Nx*Ny + (j-6)*Nx + i

#define GET3DINDEX_XYZ_Zp1 (k+1)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp2 (k+2)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp3 (k+3)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp4 (k+4)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp5 (k+5)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zp6 (k+6)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm1 (k-1)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm2 (k-2)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm3 (k-3)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm4 (k-4)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm5 (k-5)*Nx*Ny + j*Nx + i
#define GET3DINDEX_XYZ_Zm6 (k-6)*Nx*Ny + j*Nx + i

//Loops
#define FOR_X for(int i = 0; i < Nx; i++)
#define FOR_Y for(int j = 0; j < Ny; j++)
#define FOR_Z for(int k = 0; k < Nz; k++)
#define FOR_XYZ for(int ip = 0; ip < Nx*Ny*Nz; ip++) 


//Need to close off each of these macros with their corresponding END_FOR macros'
#define FOR_X0 for(int k = 0; k < Nz; k++){for(int j = 0; j < Ny; j++){int i = 0; int ip = GET3DINDEX_XYZ; 
#define FOR_X1 for(int k = 0; k < Nz; k++){for(int j = 0; j < Ny; j++){int i = Nx-1; int ip = GET3DINDEX_XYZ; 
#define FOR_Y0 for(int k = 0; k < Nz; k++){for(int i = 0; i < Nx; i++){int j = 0; int ip = GET3DINDEX_XYZ; 
#define FOR_Y1 for(int k = 0; k < Nz; k++){for(int i = 0; i < Nx; i++){int j = Ny-1; int ip = GET3DINDEX_XYZ; 
#define FOR_Z0 for(int j = 0; j < Ny; j++){for(int i = 0; i < Nx; i++){int k = 0; int ip = GET3DINDEX_XYZ; 
#define FOR_Z1 for(int j = 0; j < Ny; j++){for(int i = 0; i < Nx; i++){int k = Nz-1; int ip = GET3DINDEX_XYZ; 

#define END_FORX0 }}
#define END_FORX1 }}
#define END_FORY0 }}
#define END_FORY1 }}
#define END_FORZ0 }}
#define END_FORZ1 }}

#endif
