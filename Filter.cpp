#include "Filter.hpp"

using namespace std;

void Filter::calcFilterCoefficients(double *alphaFdir){

    //Assuming constant coefficients right now
    
    int i = 0;
    //Node 1
    b0 =      (255.0  +      alphaFdir[i])     /256.0;
    b1 =      (1.0    + 31.0*alphaFdir[i])     /32.0;
    b2 = 7.0* (-1.0   +      alphaFdir[i])     /64.0;
    b3 = 7.0* (1.0    -      alphaFdir[i])     /32.0;
    b4 = 35.0*(              alphaFdir[i]-1.0) /128.0;
    b5 = 7.0* (1.0    -      alphaFdir[i])     /32.0;
    b6 = 7.0* (-1.0   +      alphaFdir[i])     /64.0;
    b7 =      (1.0    -      alphaFdir[i])     /32.0;
    b8 =      (              alphaFdir[i]-1.0) /256.0;

    i = 1;
    //Node 2
    c0 = 127.0*alphaFdir[i]/128.0 + 1.0/256.0;
    c1 = (31.0+2.0*alphaFdir[i])/32.0;
    c2 = (7.0+50.0*alphaFdir[i])/64.0;
    c3 = 7.0*(2.0*alphaFdir[i]-1.0)/32.0;
    c4 = 35.0*(1.0-2.0*alphaFdir[i])/128.0;
    c5 = 7.0*(2.0*alphaFdir[i]-1.0)/32.0;
    c6 = 7.0*(1.0-2.0*alphaFdir[i])/64.0;
    c7 = (2.0*alphaFdir[i]-1.0)/32.0;
    c8 = (1.0-2.0*alphaFdir[i])/256.0;

    i = 2;
    //Node 3
    d0 = (2.0*alphaFdir[i] -1.0)/256.0;
    d1 = (1.0+30.0*alphaFdir[i])/32.0;
    d2 = (57.0+14.0*alphaFdir[i])/64.0;
    d3 = (7.0+18.0*alphaFdir[i])/32.0;
    d4 = 35.0*(2.0*alphaFdir[i]-1.0)/128.0;
    d5 = 7.0*(1.0-2.0*alphaFdir[i])/32.0;
    d6 = 7.0*(2.0*alphaFdir[i]-1.0)/64.0;
    d7 = (1.0-2.0*alphaFdir[i])/32.0;
    d8 = (2.0*alphaFdir[i]-1.0)/256.0;

    i = 3;
    //Node 4
    e0 = (1.0-2.0*alphaFdir[i])/256.0;
    e1 = (2.0*alphaFdir[i]-1.0)/32.0;
    e2 = (7.0+50.0*alphaFdir[i])/64.0;
    e3 = (25.0+14.0*alphaFdir[i])/32.0;
    e4 = (35.0+58.0*alphaFdir[i])/128.0;
    e5 = 7.0*(2.0*alphaFdir[i]-1.0)/32.0;
    e6 = 7.0*(1.0-2.0*alphaFdir[i])/64.0;
    e7 = (2.0*alphaFdir[i]-1.0)/32.0;
    e8 = (1.0-2.0*alphaFdir[i])/256.0;

}


void Filter::multRHSFilter(double *phi, int N, double *RHSvec){

    double cc0 = a0;
    double cc1 = a1/2.0;
    double cc2 = a2/2.0;
    double cc3 = a3/2.0;
    double cc4 = a4/2.0;
    double cc5 = a5/2.0;

    for(int ip = 0; ip < N; ip++){
        if(ip == 0){
            RHSvec[ip] = cc5*phi[N-5]  + cc4*phi[N-4]  + cc3*phi[N-3]  +
                         cc2*phi[N-2]  + cc1*phi[N-1]  + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 1){
            RHSvec[ip] = cc5*phi[N-4]  + cc4*phi[N-3]  + cc3*phi[N-2]  +
                         cc2*phi[N-1]  + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 2){
            RHSvec[ip] = cc5*phi[N-3]  + cc4*phi[N-2]  + cc3*phi[N-1]  +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 3){
            RHSvec[ip] = cc5*phi[N-2]  + cc4*phi[N-1]  + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == 4){
            RHSvec[ip] = cc5*phi[N-1]  + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }else if(ip == N-5){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[0];
        }else if(ip == N-5){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[0];
        }else if(ip == N-4){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[0]    + cc5*phi[1];
        }else if(ip == N-3){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[0]    +
                         cc4*phi[1]    + cc5*phi[2];
        }else if(ip == N-2){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[0]    + cc3*phi[1]    +
                         cc4*phi[2]    + cc5*phi[3];   
        }else if(ip == N-1){
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[0]    + cc2*phi[1]    + cc3*phi[2]    +
                         cc4*phi[3]    + cc5*phi[4];
        }else{           
            RHSvec[ip] = cc5*phi[ip-5] + cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4] + cc5*phi[ip+5];
        }                
        
    }

    
}

void Filter::FilterPeriodicY(double *phi, double *phiF){

    double RHSvec[Ny];
//    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
        double *phiPointer = &phi[ip*Ny];
        multRHSFilter(phiPointer, Ny, RHSvec);
        cyclic(offFy1, diagFy, offFy2, alphaF, alphaF, RHSvec, Ny, &phiF[ip*Ny]);
    }
}

void Filter::FilterPeriodicX(double *phi, double *phiF){

    double *phiTrans = new double[Nx*Ny];
    double *phiFTrans = new double[Nx*Ny];
    transposeMatrix(phi, Nx, Ny, phiTrans);
    //transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);

    double RHSvec[Nx];
    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Ny; ip++){
        double *phiPointer = &phiTrans[ip*Nx];
        multRHSFilter(phiPointer, Nx, RHSvec);
        cyclic(offFx1, diagFx, offFx2, alphaF, alphaF, RHSvec, Nx, &phiFTrans[ip*Nx]);
    }
    delete[] phiTrans;

    transposeMatrix(phiFTrans, Ny, Nx, phiF);
    //transposeMatrix_Fast2(phiFTrans, Ny, Nx, phiF, 60);
    delete[] phiFTrans;
}


void Filter::multRHSFilterFiniteDomain(double *phi, int N, double *RHSvec){

    double cc0 = a0;
    double cc1 = a1/2.0;
    double cc2 = a2/2.0;
    double cc3 = a3/2.0;
    double cc4 = a4/2.0;

    RHSvec[0] = b0*phi[0] + b1*phi[1] + b2*phi[2] + b3*phi[3] +
		b4*phi[4] + b5*phi[5] + b6*phi[6] + b7*phi[7] + b8*phi[8];

    RHSvec[1] = c0*phi[0] + c1*phi[1] + c2*phi[2] + c3*phi[3] +
		c4*phi[4] + c5*phi[5] + c6*phi[6] + c7*phi[7] + c8*phi[8];

    RHSvec[2] = d0*phi[0] + d1*phi[1] + d2*phi[2] + d3*phi[3] +
		d4*phi[4] + d5*phi[5] + d6*phi[6] + d7*phi[7] + d8*phi[8];

    RHSvec[3] = e0*phi[0] + e1*phi[1] + e2*phi[2] + e3*phi[3] +
		e4*phi[4] + e5*phi[5] + e6*phi[6] + e7*phi[7] + e8*phi[8];

    for(int ip = 4; ip < N-4; ip++){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
    }

    RHSvec[N-4] = e0*phi[N-1] + e1*phi[N-2] + e2*phi[N-3] + e3*phi[N-4] +
		  e4*phi[N-5] + e5*phi[N-6] + e6*phi[N-7] + e7*phi[N-8] + e8*phi[N-9];


    RHSvec[N-3] = d0*phi[N-1] + d1*phi[N-2] + d2*phi[N-3] + d3*phi[N-4] +
		  d4*phi[N-5] + d5*phi[N-6] + d6*phi[N-7] + d7*phi[N-8] + d8*phi[N-9];

    RHSvec[N-2] = c0*phi[N-1] + c1*phi[N-2] + c2*phi[N-3] + c3*phi[N-4] +
		  c4*phi[N-5] + c5*phi[N-6] + c6*phi[N-7] + c7*phi[N-8] + c8*phi[N-9];

    RHSvec[N-1] = b0*phi[N-1] + b1*phi[N-2] + b2*phi[N-3] + b3*phi[N-4] +
		  b4*phi[N-5] + b5*phi[N-6] + b6*phi[N-7] + b7*phi[N-8] + b8*phi[N-9];

};


void Filter::FilterFiniteDomainY(double *phi, double *phiF){
    
    calcFilterCoefficients(alphaFy);

    double RHSvec[Ny];
    double *work = new double[Ny];

//    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Nx; ip++){
        double *phiPointer = &phi[ip*Ny];
        multRHSFilterFiniteDomain(phiPointer, Ny, RHSvec);
        solveTri(offFy1, diagFy, offFy2, RHSvec, &phiF[ip*Ny], work, Ny);
    }

    delete[] work;

};

void Filter::FilterFiniteDomainX(double *phi, double *phiF){


    calcFilterCoefficients(alphaFx);

    double *phiTrans = new double[Nx*Ny];

    double *phiFTrans = new double[Nx*Ny];
    transposeMatrix(phi, Nx, Ny, phiTrans);
    //transposeMatrix_Fast2(phi, Nx, Ny, phiTrans, 60);

    double RHSvec[Nx];
    double *work = new double[Nx];

    #pragma omp parallel for private(RHSvec)
    for(int ip = 0 ; ip < Ny; ip++){
        double *phiPointer = &phiTrans[ip*Nx];
        multRHSFilterFiniteDomain(phiPointer, Nx, RHSvec);
        solveTri(offFx1, diagFx, offFx2, RHSvec, &phiFTrans[ip*Nx], work, Nx);
    }
/*
    for(int jp = 0; jp < Nx; jp++){
       cout << phiTrans[jp] << " " << phi[jp*Ny] << " " << endl;
    }
*/	

    delete[] phiTrans;
    delete[] work;

    transposeMatrix(phiFTrans, Ny, Nx, phiF);
    //transposeMatrix_Fast2(phiFTrans, Ny, Nx, phiF, 60);
    delete[] phiFTrans;

};
