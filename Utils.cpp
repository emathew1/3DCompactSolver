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

    #pragma omp parallel for collapse(3) num_threads(2) 
    for(int ip = 0; ip < Nx; ip++){
	for(int kp = 0; kp < Nz; kp++){
	    for(int jp = 0; jp < Ny; jp++){
		out[ip*Nz*Ny + kp*Ny + jp] = in[kp*Ny*Nx + jp*Nx + ip];
	    }
	}
    }

} 

void transposeYZXtoZXY(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for collapse(3) num_threads(2)
    for(int jp = 0; jp < Ny; jp++){
    	for(int ip = 0; ip < Nx; ip++){
	    for(int kp = 0; kp < Nz; kp++){
		out[jp*Nx*Nz + ip*Nz + kp] = in[ip*Nz*Ny + kp*Ny + jp];
	    }
	}
    }

}

void transposeXYZtoZXY(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for collapse(3) num_threads(2)
    for(int jp = 0; jp < Ny; jp++){
    	for(int ip = 0; ip < Nx; ip++){
	    for(int kp = 0; kp < Nz; kp++){
		out[jp*Nx*Nz + ip*Nz + kp] = in[kp*Nx*Ny + jp*Nx + ip];
	    }
	}
    }

}


void transposeZXYtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for collapse(3) num_threads(2) 
    for(int kp = 0; kp < Nz; kp++){
        for(int jp = 0; jp < Ny; jp++){
    	    for(int ip = 0; ip < Nx; ip++){
		out[kp*Nx*Ny + jp*Nx + ip] = in[jp*Nx*Nz + ip*Nz + kp];
	    }
	}
    }

}

void transposeYZXtoXYZ(const double *in, int Nx, int Ny, int Nz, double *out){

    #pragma omp parallel for collapse(3) num_threads(2)
    for(int kp = 0; kp < Nz; kp++){
    	for(int jp = 0; jp < Ny; jp++){
	    for(int ip = 0; ip < Nx; ip++){
		out[kp*Nx*Ny + jp*Nx + ip] = in[ip*Nz*Ny + kp*Ny + jp];
	    }
	}
    }
}


void getRange(double *phi, std::string dataName, int Nx, int Ny){
    double dataMin = 1000000;
    double dataMax = -1000000;
    for(int ip = 0; ip < Nx*Ny; ip++){
        if(phi[ip] > dataMax){
            dataMax = phi[ip];
        }

        if(phi[ip] < dataMin){
            dataMin = phi[ip];
        }
    }

    cout << "  Range of " << dataName << ": " << dataMin << ":" << dataMax << endl;

}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
