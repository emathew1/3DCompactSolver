#include <omp.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

#include "Macros.hpp"
#include "Utils.hpp"
#include "Domain.hpp"
#include "Derivatives.hpp"
#include "BC.hpp"
#include "IdealGas.hpp"
#include "VisitWriter.hpp"

int main(int argc, char *argv[]){


    cout << endl;
    cout << "----------------------------------" << endl;
    cout << " 3D 6th-Order Compressible Solver " << endl;
    cout << "     Post-Processing Edition      " << endl;
    cout << "----------------------------------" << endl;
    cout << endl;

    #pragma omp parallel
    {
      int threadCount = omp_get_num_threads();
      int threadID    = omp_get_thread_num();
      if(threadCount > 1 && threadID == 0){
          cout << " > Running with OpenMP using " << threadCount << " threads" << endl;
          cout << " > Thread limit: " << omp_get_thread_limit() << endl;
      }
    }
    omp_set_max_active_levels(2);


    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 512,
           Ny = 256,
           Nz = 128;
    double Lx = 345.0,
           Ly = 172.0,
           Lz = 86.0;;
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);

    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////
    BC::BCType bcXType = BC::PERIODIC_SOLVE;
    BC::BCType bcYType = BC::DIRICHLET_SOLVE;
    BC::BCType bcZType = BC::PERIODIC_SOLVE;

    BC::BCKind bcX0 = BC::PERIODIC;
    BC::BCKind bcX1 = BC::PERIODIC;
    BC::BCKind bcY0 = BC::SPONGE;
    BC::BCKind bcY1 = BC::SPONGE;
    BC::BCKind bcZ0 = BC::PERIODIC;
    BC::BCKind bcZ1 = BC::PERIODIC;

    BC *bc = new BC(bcXType, bcX0, bcX1,
                    bcYType, bcY0, bcY1,
                    bcZType, bcZ0, bcZ1);

    
    //create the ideal gas object...
    double mu_ref = 0.0035;
    IdealGas *ig = new IdealGas(dom, mu_ref);

    //Initialize the derivative objects we need...
    Derivatives *derivX, *derivY, *derivZ;

    cout << "  > Initializing the derivative objects " << endl;
    derivX = new Derivatives(dom, bc->bcXType, Derivatives::DIRX);
    derivY = new Derivatives(dom, bc->bcYType, Derivatives::DIRY);
    derivZ = new Derivatives(dom, bc->bcZType, Derivatives::DIRZ);

    //Initialize the data containers that we'll need...
    double *rho, *rhoU, *rhoV, *rhoW, *rhoE;
    double *U, *V, *W, *T, *p;

    cout << "  > Allocating base memory " << endl;
    rho  = new double[Nx*Ny*Nz];
    rhoU = new double[Nx*Ny*Nz];
    rhoV = new double[Nx*Ny*Nz];
    rhoW = new double[Nx*Ny*Nz];
    rhoE = new double[Nx*Ny*Nz];

    U = new double[Nx*Ny*Nz];
    V = new double[Nx*Ny*Nz];
    W = new double[Nx*Ny*Nz];
    T = new double[Nx*Ny*Nz];
    p = new double[Nx*Ny*Nz];

    
    double temp;

    //loading density file
    cout << "  > Loading density file...";
    ifstream infile;
    string densityFile = "ShearLayer/A3_FullDomain_k6nocurl_10oversqrt3percpert/rho.out.5000";
    infile.open(densityFile);
    if(!infile){
	cout << "Unable to open density file: " << densityFile << endl;
	exit(1); 
    }    
    int count = 0;
    while(infile >> temp){
	rho[count] = temp;
	count++;
    }
    cout << "done!" << endl;
    infile.close();
    getRange(rho, "RHO", Nx, Ny, Nz);

    //loading x momentum file
    cout << "  > Loading rhoU file...";
    string rhoUFile = "ShearLayer/A3_FullDomain_k6nocurl_10oversqrt3percpert/rhoU.out.5000";
    infile.open(rhoUFile);
    if(!infile){
	cout << "Unable to open rhoU file: " << rhoUFile << endl;
	exit(1); 
    }    
    count = 0;
    while(infile >> temp){
	rhoU[count] = temp;
	count++;
    }
    cout << "done!" << endl;
    infile.close();
    getRange(rhoU, "RHOU", Nx, Ny, Nz);

    //loading y momentum file
    cout << "  > Loading rhoV file...";
    string rhoVFile = "ShearLayer/A3_FullDomain_k6nocurl_10oversqrt3percpert/rhoV.out.5000";
    infile.open(rhoVFile);
    if(!infile){
	cout << "Unable to open rhoV file: " << rhoVFile << endl;
	exit(1); 
    }    
    count = 0;
    while(infile >> temp){
	rhoV[count] = temp;
	count++;
    }
    cout << "done!" << endl;
    infile.close();

    getRange(rhoV, "RHOV", Nx, Ny, Nz);

    //loading z momentum file
    cout << "  > Loading rhoW file...";
    string rhoWFile = "ShearLayer/A3_FullDomain_k6nocurl_10oversqrt3percpert/rhoW.out.5000";
    infile.open(rhoWFile);
    if(!infile){
	cout << "Unable to open rhoW file: " << rhoWFile << endl;
	exit(1); 
    }    
    count = 0;
    while(infile >> temp){
	rhoW[count] = temp;
	count++;
    }
    cout << "done!" << endl;
    infile.close();

    //Loading energy file
    cout << "  > Loading rhoE file...";
    string rhoEFile = "ShearLayer/A3_FullDomain_k6nocurl_10oversqrt3percpert/rhoE.out.5000";
    infile.open(rhoEFile);
    if(!infile){
	cout << "Unable to open rhoE file: " << rhoEFile << endl;
	exit(1); 
    }    
    count = 0;
    while(infile >> temp){
	rhoE[count] = temp;
	count++;
    }
    cout << "done!" << endl;
    infile.close();
    getRange(rhoE, "RHOE", Nx, Ny, Nz);


    cout << "  > Calculating ideal gas data from the input files." << endl;

    #pragma omp parallel for
    FOR_XYZ{
	U[ip] = rhoU[ip]/rho[ip];
	V[ip] = rhoV[ip]/rho[ip];
	W[ip] = rhoW[ip]/rho[ip];
	p[ip] = ig->solvep(rho[ip], rhoE[ip], U[ip], V[ip], W[ip]); 
	T[ip] = ig->solveT(rho[ip], p[ip]);		
    }
    cout << "  > Done!" << endl;

    //Calculating cross derivatives of velocity to get vorticity...
    double *Uy, *Uz;
    double *Vx, *Vz;
    double *Wx, *Wy;
    double *transTemp1, *transTemp2;

    cout << "  > Allocating memory to calculate vorticity..." << endl;
    Uy = new double[Nx*Ny*Nz];
    Uz = new double[Nx*Ny*Nz];
    Vx = new double[Nx*Ny*Nz];
    Vz = new double[Nx*Ny*Nz];
    Wx = new double[Nx*Ny*Nz];
    Wy = new double[Nx*Ny*Nz];
    transTemp1 = new double[Nx*Ny*Nz];
    transTemp2 = new double[Nx*Ny*Nz];
    cout << "  > Done!" << endl;

    int blocksize = 16;
    
    //X-Derivatives...
    cout << "  > X-Derivatives..." << endl;
    derivX->calc1stDerivField(V, Vx);
    derivX->calc1stDerivField(W, Wx);

    
    cout << "  > Y-Derivatives..." << endl;
    transposeXYZtoYZX_Fast(U, Nx, Ny, Nz, transTemp1, blocksize);
    derivY->calc1stDerivField(transTemp1, transTemp2);
    transposeYZXtoXYZ_Fast(transTemp2, Nx, Ny, Nz, Uy, blocksize);

    transposeXYZtoYZX_Fast(W, Nx, Ny, Nz, transTemp1, blocksize);
    derivY->calc1stDerivField(transTemp1, transTemp2);
    transposeYZXtoXYZ_Fast(transTemp2, Nx, Ny, Nz, Wy, blocksize);

    cout << "  > Z-Derivatives..." << endl;
    transposeXYZtoZXY_Fast(U, Nx, Ny, Nz, transTemp1, blocksize);
    derivZ->calc1stDerivField(transTemp1, transTemp2);
    transposeZXYtoXYZ_Fast(transTemp2, Nx, Ny, Nz, Uz, blocksize);

    transposeXYZtoZXY_Fast(V, Nx, Ny, Nz, transTemp1, blocksize);
    derivZ->calc1stDerivField(transTemp1, transTemp2);
    transposeZXYtoXYZ_Fast(transTemp2, Nx, Ny, Nz, Vz, blocksize);
    cout << " > Done! " << endl;

    delete[] transTemp1;
    delete[] transTemp2;

    double *vortMag = new double[Nx*Ny*Nz];
    
    #pragma omp parallel for
    FOR_XYZ{
	double vortx = Wy[ip] - Vz[ip];
	double vorty = Uz[ip] - Wx[ip];
	double vortz = Vx[ip] - Uy[ip];
	vortMag[ip] = sqrt(vortx*vortx + vorty*vorty + vortz*vortz);
    }

    getRange(vortMag, "VORTMAG", Nx, Ny, Nz);


    //output needs to be in float format...
    float *rhof = new float[Nx*Ny*Nz];
    float *pf = new float[Nx*Ny*Nz];
    float *Tf = new float[Nx*Ny*Nz];
    float *Uf = new float[Nx*Ny*Nz];
    float *Vf = new float[Nx*Ny*Nz];
    float *Wf = new float[Nx*Ny*Nz];
    float *vortMagf = new float[Nx*Ny*Nz];
  
    #pragma omp parallel for
    FOR_XYZ{
	rhof[ip] = (float)rho[ip];
	pf[ip]   = (float)p[ip];
	Tf[ip]   = (float)T[ip];
	Uf[ip]   = (float)U[ip];
	Vf[ip]   = (float)V[ip];
	Wf[ip] 	 = (float)W[ip];
	vortMagf[ip] = (float)vortMag[ip];

    }
   

    const int numvars = 7;
    int dims[] = {Nx, Ny, Nz};
    int vardims[] = {1, 1, 1, 1, 1, 1, 1};
    int centering[] = {1, 1, 1, 1, 1, 1, 1};
    const char * const varnames[] = {"RHO", "P", "U", "V", "W", "T", "VORTMAG"};
    float *arrays[] = {(float*)rhof, (float*)pf, (float*)Uf, (float*)Vf, (float*)Wf, (float*)Tf, (float*)vortMagf};
    write_regular_mesh("test.vtk", !0, dims, numvars, vardims,  centering, varnames, arrays);

    return 0;
}
