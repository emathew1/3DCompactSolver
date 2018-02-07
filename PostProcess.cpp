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

    
    //Initialize the derivative objects we need...
    Derivatives *derivX, *derivY, *derivZ;

    derivX = new Derivatives(dom, bc->bcXType, Derivatives::DIRX);
    derivY = new Derivatives(dom, bc->bcYType, Derivatives::DIRY);
    derivZ = new Derivatives(dom, bc->bcZType, Derivatives::DIRZ);

    //Initialize the data containers that we'll need...
    double *rho, *rhoU, *rhoV, *rhoW, *rhoE;
    double *U, *V, *W, *T, *p;

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
    cout << "Loading density file...";
    ifstream infile;
    string densityFile = "rho.out.5000";
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

    //loading x momentum file
    cout << "Loading rhoU file...";
    string rhoUFile = "rhoU.out.5000";
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

    //loading y momentum file
    cout << "Loading rhoV file...";
    string rhoVFile = "rhoV.out.5000";
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

    //loading z momentum file
    cout << "Loading rhoW file...";
    string rhoWFile = "rhoW.out.5000";
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
    cout << "Loading rhoE file...";
    string rhoEFile = "rhoE.out.5000";
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



    return 0;
}
