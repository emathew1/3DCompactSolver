#include "Filter.hpp"

using namespace std;

void Filter::multRHSPeriodicFilter(double *phi, double *RHSvec){

    double cc0 = a0_8;
    double cc1 = a1_8/2.0;
    double cc2 = a2_8/2.0;
    double cc3 = a3_8/2.0;
    double cc4 = a4_8/2.0;

    for(int ip = 0; ip < N; ip++){
        if(ip == 0){
            RHSvec[ip] = cc4*phi[N-4]  + cc3*phi[N-3]  +
                         cc2*phi[N-2]  + cc1*phi[N-1]  + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == 1){
            RHSvec[ip] = cc4*phi[N-3]  + cc3*phi[N-2]  +
                         cc2*phi[N-1]  + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == 2){
            RHSvec[ip] = cc4*phi[N-2]  + cc3*phi[N-1]  +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == 3){
            RHSvec[ip] = cc4*phi[N-1]  + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }else if(ip == N-4){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[0];
        }else if(ip == N-3){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[0]    +
                         cc4*phi[1];
        }else if(ip == N-2){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[0]    + cc3*phi[1]    +
                         cc4*phi[2];   
        }else if(ip == N-1){
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[0]    + cc2*phi[1]    + cc3*phi[2]    +
                         cc4*phi[3];
        }else{           
            RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                         cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                         cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                         cc4*phi[ip+4];
        }                
        
    }
}

void Filter::multRHSDirichletFilter(double *phi, double *RHSvec){

    double aa0 = a0_6;
    double aa1 = a1_6/2.0;
    double aa2 = a1_6/2.0;
    double aa3 = a1_6/2.0;

    double cc0 = a0_8;
    double cc1 = a1_8/2.0;
    double cc2 = a2_8/2.0;
    double cc3 = a3_8/2.0;
    double cc4 = a4_8/2.0;

    RHSvec[0] = a00*phi[0] + a01*phi[1] + a02*phi[2] + a03*phi[3] + a04*phi[4];
    RHSvec[1] = b00*phi[0] + b01*phi[1] + b02*phi[2] + b03*phi[3] + b04*phi[4];
    RHSvec[2] = c00*phi[0] + c01*phi[1] + c02*phi[2] + c03*phi[3] + c04*phi[4];
    RHSvec[3] = aa3*phi[0] + aa2*phi[1] + aa1*phi[2] + aa0*phi[3] + aa1*phi[4] + 
		aa2*phi[5] + aa3*phi[6];
    
    RHSvec[N-1] = a00*phi[N-1] + a01*phi[N-2] + a02*phi[N-3] + a03*phi[N-4] + a04*phi[N-5];
    RHSvec[N-2] = b00*phi[N-1] + b01*phi[N-2] + b02*phi[N-3] + b03*phi[N-4] + b04*phi[N-5];
    RHSvec[N-3] = c00*phi[N-1] + c01*phi[N-2] + c02*phi[N-3] + c03*phi[N-4] + c04*phi[N-5];
    RHSvec[N-4] = aa3*phi[N-1] + aa2*phi[N-2] + aa1*phi[N-3] + aa0*phi[N-4] + aa1*phi[N-5] + 
		  aa2*phi[N-6] + aa3*phi[N-7];
 

    for(int ip = 4; ip < N-4; ip++){
        RHSvec[ip] = cc4*phi[ip-4] + cc3*phi[ip-3] +
                     cc2*phi[ip-2] + cc1*phi[ip-1] + cc0*phi[ip]   +
                     cc1*phi[ip+1] + cc2*phi[ip+2] + cc3*phi[ip+3] +
                     cc4*phi[ip+4];
    }
}

void Filter::FilterPeriodic(double *phi, double *phiF){

    double RHSvec[Ny];
    multRHSPeriodicFilter(phi, RHSvec);
    cyclic(offlowerF, diagF, offupperF, alphaF, alphaF, RHSvec, N, phiF);

}

void Filter::FilterDirichlet(double *phi, double *phiF){
    

    double RHSvec[Ny];
    double *work = new double[Ny];

    multRHSDirichletFilter(phi, RHSvec);
    solveTri(offlowerF, diagF, offupperF, RHSvec, phiF, work, N);

    delete[] work;

}

void Filter::compactFilter(double *phi, double *phiF){

    if(bcType == BC::PERIODIC_SOLVE){
	FilterPeriodic(phi, phiF);
    }else if(bcType == BC::DIRICHLET_SOLVE){
	FilterDirichlet(phi, phiF);
    }

}
