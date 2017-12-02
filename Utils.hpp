#ifndef _UTILSH_
#define _UTILSH_

#include <math.h>
#include <cstring>
#include <iostream>

using namespace std;

void solveTri(double a[], double b[], double c[], double d[], double x[], double *work, int n);
void cyclic(double *a, double *b, double *c, double alpha, double beta,
        double *r, int n, double *x);
void transposeMatrix(double *in, int Nx, int Ny, double *out);
void transposeMatrix_Fast1(const double *in, int n, int p, double *out, int block);
void transposeMatrix_Fast2(const double *in, int n, int p, double *out, int blocksize);

void transposeXYZtoYZX(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeXYZtoYZX_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeXYZtoZXY(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeXYZtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeYZXtoZXY(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeYZXtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeZXYtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeZXYtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);
void transposeYZXtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out);
void transposeYZXtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize);


void getRange(double *phi, std::string dataName, int Nx, int Ny);
double fRand(double fMin, double fMax);

#endif
