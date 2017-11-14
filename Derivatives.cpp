#include "Derivatives.hpp"

using namespace std;

void Derivatives::multRHSDerivPeriodic(double dh, double *phi, int N, double *RHSvec){

    double c1 = -c/6.0;
    double c2 = -b/4.0;
    double c3 = -a/2.0;
    double c4 =  a/2.0;
    double c5 =  b/4.0;
    double c6 =  c/6.0;

    RHSvec[0] = c1*phi[N-3] + c2*phi[N-2] + c3*phi[N-1] + \
                        c4*phi[1] + c5*phi[2] + c6*phi[3];
    RHSvec[1] = c1*phi[N-2] + c2*phi[N-1] + c3*phi[0] + \
                        c4*phi[2] + c5*phi[3] + c6*phi[4];
    RHSvec[2] = c1*phi[N-1] + c2*phi[0] + c3*phi[1] + \
                        c4*phi[3] + c5*phi[4] + c6*phi[5];

    for(int ip = 3; ip < N-3; ip++){
        RHSvec[ip] = c1*phi[ip-3] + c2*phi[ip-2] + c3*phi[ip-1] + \
                        c4*phi[ip+1] + c5*phi[ip+2] + c6*phi[ip+3];

    }

    RHSvec[N-3] = c1*phi[N-6] + c2*phi[N-5] + c3*phi[N-4] + \
                        c4*phi[N-2] + c5*phi[N-1] + c6*phi[0];
    RHSvec[N-2] = c1*phi[N-5] + c2*phi[N-4] + c3*phi[N-3] + \
                        c4*phi[N-1] + c5*phi[0] + c6*phi[1];
    RHSvec[N-1] = c1*phi[N-4] + c2*phi[N-3] + c3*phi[N-2] + \
                        c4*phi[0] + c5*phi[1] + c6*phi[2];

    for(int ip = 0; ip < N; ip++){
        RHSvec[ip] /= dh;
    }

}

void Derivatives::CompactDYPeriodic(double *phi, double *dphidy){

    double RHSvec[Ny];
//    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
        double *phiPointer = &phi[ip*Ny];
        multRHSDerivPeriodic(dy, phiPointer, Ny, RHSvec);
        cyclic(offy1, diagy, offy2, alpha, alpha, RHSvec, Ny, &dphidy[ip*Ny]);
    }

}

void Derivatives::CompactDXPeriodic(double *phi, double *dphidx){

    double *phiTrans = new double[Nx*Ny];

    double *dphidxTrans = new double[Nx*Ny];
    transposeMatrix(phi, Nx, Ny, phiTrans);
    //transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);


    double RHSvec[Nx];
    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Ny; ip++){
        double *phiPointer = &phiTrans[ip*Nx];
        multRHSDerivPeriodic(dx, phiPointer, Nx, RHSvec);
        cyclic(offx1, diagx, offx2, alpha, alpha, RHSvec, Nx, &dphidxTrans[ip*Nx]);
    }
    delete[] phiTrans;

    transposeMatrix(dphidxTrans, Ny, Nx, dphidx);
    //transposeMatrix_Fast2(dphidxTrans, Ny, Nx, dphidx, 60);
    delete[] dphidxTrans;
}

void Derivatives::multRHSDerivDirichlet(double dh, double *phi, int N, double *RHSvec){


    double cc1 = -c/6.0;
    double cc2 = -b/4.0;
    double cc3 = -a/2.0;
    double cc4 =  a/2.0;
    double cc5 =  b/4.0;
    double cc6 =  c/6.0;

    RHSvec[0] = a1*phi[0] + b1*phi[1] + c1*phi[2] + \
                        d1*phi[3] + e1*phi[4] + f1*phi[5];
    RHSvec[1] = a2*phi[0] + b2*phi[1] + c2*phi[2] + \
                        d2*phi[3] + e2*phi[4];

    for(int ip = 2; ip < N-2; ip++){
        RHSvec[ip] = cc2*phi[ip-2] + cc3*phi[ip-1] + \
                        cc4*phi[ip+1] + cc5*phi[ip+2];

    }

    RHSvec[N-2] = -e2*phi[N-5] - d2*phi[N-4] - \
                        c2*phi[N-3] - b2*phi[N-2] - a2*phi[N-1];
    RHSvec[N-1] = -f1*phi[N-6] - e1*phi[N-5] - d1*phi[N-4] - \
                        c1*phi[N-3] - b1*phi[N-2] - a1*phi[N-1];

    for(int ip = 0; ip < N; ip++){
        RHSvec[ip] /= dh;
    }

}

void Derivatives::CompactDYDirichlet(double *phi, double *dphidy){

    double RHSvec[Ny];
    double *work = new double[Ny];

//    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
        double *phiPointer = &phi[ip*Ny];
        multRHSDerivDirichlet(dy, phiPointer, Ny, RHSvec);
        solveTri(offy1, diagy, offy2, RHSvec, &dphidy[ip*Ny], work, Ny);
    }

    delete[] work;
}

void Derivatives::CompactDXDirichlet(double *phi, double *dphidx){

    double *phiTrans = new double[Nx*Ny];

    double *dphidxTrans = new double[Nx*Ny];
    transposeMatrix(phi, Nx, Ny, phiTrans);
    //transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);


    double RHSvec[Nx];
    double *work = new double[Nx];

    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Ny; ip++){
        double *phiPointer = &phiTrans[ip*Nx];
        multRHSDerivDirichlet(dx, phiPointer, Nx, RHSvec);
        solveTri(offx1, diagx, offx2, RHSvec, &dphidxTrans[ip*Nx], work, Nx);
    }
    delete[] phiTrans;
    delete[] work;

    transposeMatrix(dphidxTrans, Ny, Nx, dphidx);
    //transposeMatrix_Fast2(dphidxTrans, Ny, Nx, dphidx, 60);
    delete[] dphidxTrans;
}

