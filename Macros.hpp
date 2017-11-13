#ifndef MACROS_HPP
#define MACROS_HPP

//Need to be sure that Nx,Ny,Nz are defined within scope to use these
#define FOR_X for(int i = 0; i < Nx; i++)
#define FOR_Y for(int j = 0; j < Ny; j++)
#define FOR_Z for(int k = 0; k < Nz; k++)
#define FOR_XYZ for(int ip = 0; ip < Nx*Ny*Nz; ip++) 

#define GET3DINDEX k*Nx*Ny + j*Nx + i

#endif
