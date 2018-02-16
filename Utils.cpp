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

void getBaseNodeIndex(Domain *dom, double xp[3], int (&ind)[3]){

    if( xp[0] < 0.0 || xp[0] > dom->Lx || 
	xp[1] < 0.0 || xp[1] > dom->Ly || 
	xp[2] < 0.0 || xp[2] > dom->Lz){

	cout << "Error: in findNearestNodeIndices point entered was out of bounds!" << endl;
	cout << "xp = { " << xp[0] << ", " << xp[1] << ", " << xp[2] << "}" << endl;
	cout << "dom = { " << dom->Lx << ", " << dom->Ly << ", " << dom->Lz << "}" << endl;
    }else{

	ind[0] = xp[0]/dom->dx;
	ind[1] = xp[1]/dom->dy;
	ind[2] = xp[2]/dom->dz;

	cout << "ind = {" << ind[0] << " " << ind[1] << " " << ind[2] <<"}" << endl;
    }

}

double linearInterpolation(Domain *dom, double *fieldIn, double xp[3]){

    int baseIndex[3];
    double Nx = dom->Nx, Ny = dom->Ny, Nz = dom->Nz;
    double c000, c001, c010, c011, c100, c101, c110, c111;
    double c00, c01, c10, c11, c0, c1;
    int    i000, i001, i010, i011, i100, i101, i110, i111;
    int    i, j, k;
    double x0, x1, y0, y1, z0, z1, xd, yd, zd;

   
    getBaseNodeIndex(dom, xp, baseIndex);

    x0 = dom->x[baseIndex[0]];
    x1 = dom->x[baseIndex[0]+1];
    y0 = dom->y[baseIndex[1]];
    y1 = dom->y[baseIndex[1]+1];
    z0 = dom->z[baseIndex[2]];
    z1 = dom->z[baseIndex[2]+1];

    xd = (xp[0] - x0)/(x1 - x0);
    yd = (xp[1] - y0)/(y1 - y0);
    zd = (xp[2] - z0)/(z1 - z0);

    i = baseIndex[0];   j = baseIndex[1];   k = baseIndex[2];    
    i000 = GET3DINDEX_XYZ;

    i = baseIndex[0];   j = baseIndex[1];   k = baseIndex[2]+1;
    i001 = GET3DINDEX_XYZ;

    i = baseIndex[0];   j = baseIndex[1]+1; k = baseIndex[2];
    i010 = GET3DINDEX_XYZ;

    i = baseIndex[0];   j = baseIndex[1]+1; k = baseIndex[2]+1;
    i011 = GET3DINDEX_XYZ;

    i = baseIndex[0]+1; j = baseIndex[1];   k = baseIndex[2];
    i100 = GET3DINDEX_XYZ;

    i = baseIndex[0]+1; j = baseIndex[1];   k = baseIndex[2]+1;
    i101 = GET3DINDEX_XYZ;

    i = baseIndex[0]+1; j = baseIndex[1]+1; k = baseIndex[2];
    i110 = GET3DINDEX_XYZ;

    i = baseIndex[0]+1; j = baseIndex[1]+1; k = baseIndex[2]+1;
    i111 = GET3DINDEX_XYZ;

    c000 = fieldIn[i000];
    c001 = fieldIn[i001];
    c010 = fieldIn[i010];
    c011 = fieldIn[i011];
    c100 = fieldIn[i100];
    c101 = fieldIn[i101];
    c110 = fieldIn[i110];
    c111 = fieldIn[i111];

    //Interpolation in the x direction first
    c00  = c000*(1 - xd) + c100*xd;
    c01  = c001*(1 - xd) + c101*xd;
    c10  = c010*(1 - xd) + c110*xd;
    c11  = c011*(1 - xd) + c111*xd;

    //Interpolation in the y direction second
    c0   = c00*(1 - yd) + c10*yd;
    c1   = c01*(1 - yd) + c11*yd;

    //Interpolation in the z direction last
    return  c0*(1 - zd) + c1*zd; 

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

void getRangeValue(double *phi, int Nx, int Ny, int Nz, double &minVal, double &maxVal){

    double dataMin = 1000000;
    #pragma omp parallel for reduction(min:dataMin)
    for(int ip = 0; ip < Nx*Ny*Nz; ip++)
        dataMin = min(dataMin, phi[ip]);

    double dataMax = -1000000;
    #pragma omp parallel for reduction(max:dataMax)
    for(int ip = 0; ip < Nx*Ny*Nz; ip++) 
        dataMax = max(dataMax, phi[ip]);
    

    minVal = dataMin;
    maxVal = dataMax;
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
