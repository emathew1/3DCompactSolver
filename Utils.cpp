#include "Utils.hpp"


void solveTri(double a[], double b[], double c[], double d[], double x[], double *work, int size)
{
        memcpy(work, b, size*sizeof(double));
        memcpy(x, d, size*sizeof(double));

        for(int ip = 1; ip < size; ip++){
            double m = a[ip]/work[ip-1];
            work[ip] = work[ip] - m*c[ip-1];
            x[ip] = x[ip] - m*x[ip-1];
        }

        x[size-1] /= work[size-1];
        for(int ip = size-2; ip >= 0; ip--){
            x[ip] = (x[ip] - c[ip]*x[ip+1])/work[ip];
        }

}

void cyclic(double *a, double *b, double *c, double alpha, double beta, double *r, int n, double *x)
{
        unsigned long i;
        double fact,gamma,bb[n],u[n],z[n],work[n];

        if (n <= 2) cout << "n too small in cyclic" << endl;

        gamma = -b[0];
        for (i=0;i<n;i++) bb[i]=b[i];
        bb[0]=b[0]-gamma;
        bb[n-1]=b[n-1]-alpha*beta/gamma;

        solveTri(a,bb,c,r,x,work,n);

        for (i=0;i<n;i++) u[i]=0.0;
        u[0]=gamma;
        u[n-1]=alpha;

        solveTri(a,bb,c,u,z,work,n);

        fact=(x[0]+beta*x[n-1]/gamma)/(1.0+z[0]+beta*z[n-1]/gamma);

        for (i=0;i<n;i++) x[i] -= fact*z[i];

}

void transposeMatrix(double *in, int Nx, int Ny, double *out){

//    #pragma omp parallel for collapse(2)
    for(int i = 0; i < Ny; ++i)
        for(int j = 0; j < Nx; ++j)
            out[i*Nx + j] =  in[j*Ny + i];

}

void transposeMatrix_Fast1(const double *in, int n, int p, double *out, int block){
//    #pragma omp parallel for
    for (int i = 0; i < n; i += block) {
        for(int j = 0; j < n; ++j) {
            for(int b = 0; b < block && i + b < n; ++b) {
                out[j*n + i + b] = in[(i + b)*n + j];
            }
        }
    }
}

void transposeMatrix_Fast2(const double *in, int n, int p, double *out, int blocksize){
    int i, j, row, col;
//    #pragma omp parallel for private(i, j, row, col) collapse(2) // schedule(static, 2)
    for ( i = 0; i < n; i += blocksize) {
        for ( j = 0; j < p; j += blocksize) {
            for (row = i; row < i + blocksize && row < n; row++) {
                for (col = j; col < j + blocksize && col < p; col++) {
                    out[row*p + col] = in[col*n + row];
                }
            }
        }
    }
}

void transposeXYZtoYZX(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for schedule(static) collapse(3) num_threads(NUMTHREADSNEST) 
    for(int ip = 0; ip < Nx; ip++){
	for(int kp = 0; kp < Nz; kp++){
	    for(int jp = 0; jp < Ny; jp++){
		out[ip*Nz*Ny + kp*Ny + jp] = in[kp*Ny*Nx + jp*Nx + ip];
	    }
	}
    }

} 

void transposeXYZtoYZX_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    #pragma omp parallel for private(i, j, k, row, col, sli) collapse(3) schedule(static) num_threads(NUMTHREADSNEST)
    for ( i = 0; i < Nx; i += blocksize) {
	for( k = 0; k < Nz; k += blocksize){
            for ( j = 0; j < Ny; j += blocksize){
                for (row = i; row < i + blocksize && row < Nx; row++) {
		    for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                    	for (col = j; col < j + blocksize && col < Ny; col++) {
                            out[row*Ny*Nz + sli*Ny + col] = in[sli*Ny*Nx + col*Nx + row];
			}
                    }
	        }
            }
        }
    }
}



void transposeYZXtoZXY(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for schedule(static) collapse(3) num_threads(NUMTHREADSNEST)
    for(int jp = 0; jp < Ny; jp++){
    	for(int ip = 0; ip < Nx; ip++){
	    for(int kp = 0; kp < Nz; kp++){
		out[jp*Nx*Nz + ip*Nz + kp] = in[ip*Nz*Ny + kp*Ny + jp];
	    }
	}
    }

}

void transposeYZXtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    #pragma omp parallel for private(i, j, k, row, col, sli) collapse(3) schedule(static) num_threads(NUMTHREADSNEST)
    for ( j = 0; j < Ny; j += blocksize){
        for ( i = 0; i < Nx; i += blocksize) {
	    for( k = 0; k < Nz; k += blocksize){
                for (col = j; col < j + blocksize && col < Ny; col++) {
                    for (row = i; row < i + blocksize && row < Nx; row++) {
		        for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                            out[col*Nx*Nz + row*Nz + sli] = in[row*Ny*Nz + sli*Ny + col];
			}
                    }
	        }
            }
        }
    }
}



void transposeXYZtoZXY(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for schedule(static) collapse(3) num_threads(NUMTHREADSNEST)
    for(int jp = 0; jp < Ny; jp++){
    	for(int ip = 0; ip < Nx; ip++){
	    for(int kp = 0; kp < Nz; kp++){
		out[jp*Nx*Nz + ip*Nz + kp] = in[kp*Nx*Ny + jp*Nx + ip];
	    }
	}
    }

}


void transposeXYZtoZXY_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    #pragma omp parallel for private(i, j, k, row, col, sli) collapse(3) schedule(static) num_threads(NUMTHREADSNEST)
    for ( j = 0; j < Ny; j += blocksize){
        for ( i = 0; i < Nx; i += blocksize) {
	    for( k = 0; k < Nz; k += blocksize){
                for (col = j; col < j + blocksize && col < Ny; col++) {
                    for (row = i; row < i + blocksize && row < Nx; row++) {
		        for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                            out[col*Nx*Nz + row*Nz + sli] = in[sli*Ny*Nx + col*Nx + row];
			}
                    }
	        }
            }
        }
    }
}


void transposeZXYtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for schedule(static) collapse(3) num_threads(NUMTHREADSNEST) 
    for(int kp = 0; kp < Nz; kp++){
        for(int jp = 0; jp < Ny; jp++){
    	    for(int ip = 0; ip < Nx; ip++){
		out[kp*Nx*Ny + jp*Nx + ip] = in[jp*Nx*Nz + ip*Nz + kp];
	    }
	}
    }

}

void transposeZXYtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    #pragma omp parallel for private(i, j, k, row, col, sli) collapse(3) schedule(static) num_threads(NUMTHREADSNEST)
    for( k = 0; k < Nz; k += blocksize){
        for ( j = 0; j < Ny; j += blocksize){
            for ( i = 0; i < Nx; i += blocksize) {
		for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                    for (col = j; col < j + blocksize && col < Ny; col++) {
                        for (row = i; row < i + blocksize && row < Nx; row++) {
                            out[sli*Nx*Ny + col*Nx + row] = in[col*Nz*Nx + row*Nz + sli];
			}
                    }
	        }
            }
        }
    }
}

void transposeYZXtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for schedule(static) collapse(3) num_threads(NUMTHREADSNEST)
    for(int kp = 0; kp < Nz; kp++){
    	for(int jp = 0; jp < Ny; jp++){
	    for(int ip = 0; ip < Nx; ip++){
		out[kp*Nx*Ny + jp*Nx + ip] = in[ip*Nz*Ny + kp*Ny + jp];
	    }
	}
    }
}

void transposeYZXtoXYZ_Fast(const double *in, int Nx, int Ny, int Nz, double *out, int blocksize){
    int i, j, k, row, col, sli;
    #pragma omp parallel for private(i, j, k, row, col, sli) collapse(3) schedule(static) num_threads(NUMTHREADSNEST)
    for( k = 0; k < Nz; k += blocksize){
        for ( j = 0; j < Ny; j += blocksize){
            for ( i = 0; i < Nx; i += blocksize) {
		for(sli = k; sli < k + blocksize && sli < Nz; sli++){
                    for (col = j; col < j + blocksize && col < Ny; col++) {
                        for (row = i; row < i + blocksize && row < Nx; row++) {
                            out[sli*Nx*Ny + col*Nx + row] = in[row*Nz*Ny + sli*Ny + col];
			}
                    }
	        }
            }
        }
    }
}



void getRange(double *phi, std::string dataName, int Nx, int Ny, int Nz){

    double dataMin = 1000000;
    #pragma omp parallel for reduction(min:dataMin)
    for(int ip = 0; ip < Nx*Ny*Nz; ip++)
        dataMin = min(dataMin, phi[ip]);

    double dataMax = -1000000;
    #pragma omp parallel for reduction(max:dataMax)
    for(int ip = 0; ip < Nx*Ny*Nz; ip++) 
        dataMax = max(dataMax, phi[ip]);
    

    cout << "  Range of " << dataName << ": " << dataMin << ":" << dataMax << endl;

}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


void tic()
{
 t0 = Clock::now();
}

void toc()
{
    Clock::time_point t1 = Clock::now();
    milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
    std::cout <<" " << std::setw(6) << ms.count() << " ms\n";
}
