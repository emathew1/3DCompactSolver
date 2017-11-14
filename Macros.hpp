#ifndef _MACROSH_
#define _MACROSH_

//Need to be sure that Nx,Ny,Nz are defined within scope to use these
#define FOR_X for(int i = 0; i < Nx; i++)
#define FOR_Y for(int j = 0; j < Ny; j++)
#define FOR_Z for(int k = 0; k < Nz; k++)
#define FOR_XYZ for(int ip = 0; ip < Nx*Ny*Nz; ip++) 

#define GET3DINDEX_XYZ k*Nx*Ny + j*Nx + i
#define GET3DINDEX_YZX i*Nz*Ny + k*Ny + j
#define GET3DINDEX_ZXY j*Nx*Nz + i*Nz + k

#endif
