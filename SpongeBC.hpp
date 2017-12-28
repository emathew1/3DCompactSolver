#ifndef _SPONGEBCH_
#define _SPONGEBCH_

#include <iostream>
#include "Macros.hpp"
#include "Domain.hpp"
#include "IdealGas.hpp"
#include "BC.hpp"

class SpongeBC{

    public:

	Domain *domain;
	IdealGas *idealGas;
	BC *bc;

	int Nx, Ny, Nz, N;

	double avgT;
	double epsP;
	double spongeP;
	double spongeStrength;
	double spongeLX;
	double spongeLY;
	double spongeLZ;

	double *sigma;
	double *spongeRhoAvg;
	double *spongeRhoUAvg;
	double *spongeRhoVAvg;
	double *spongeRhoWAvg;
	double *spongeRhoEAvg;
    
	SpongeBC(Domain *domain, IdealGas *idealGas, BC *bc){

	    std::cout << endl;
	    std::cout << " > Sponge BC found, initializing Sponge average fields and strength fields..." << std::endl;
	
	    this->domain = domain;
	    this->idealGas = idealGas;
	    this->bc = bc;

	    this->Nx = domain->Nx;
	    this->Ny = domain->Ny;
	    this->Nz = domain->Nz;
	    N = Nx*Ny*Nz;

	    sigma 	  = new double[Nx*Ny*Nz];
	    spongeRhoAvg  = new double[Nx*Ny*Nz];
	    spongeRhoUAvg = new double[Nx*Ny*Nz];
	    spongeRhoVAvg = new double[Nx*Ny*Nz];
	    spongeRhoWAvg = new double[Nx*Ny*Nz];
	    spongeRhoEAvg = new double[Nx*Ny*Nz];

	    avgT = 5.0;
	    epsP = 0.005;
	    spongeP = 1.0/idealGas->gamma;
	    spongeStrength = 12.0;
	    spongeLX = 0.25*domain->Lx;
	    spongeLY = 0.125*domain->Ly;
	    spongeLZ = 0.25*domain->Lz;

	
	    //Use this data to initialize the sponge zones / sponge sigma strength...
	    if(bc->bcX0 == BC::SPONGE){
		FOR_X{
		    if(domain->x[i] < spongeLX){
		        double spongeX = (spongeLX - domain->x[i])/spongeLX;
		        FOR_Y{
			    FOR_Z{
			        int ii = GET3DINDEX_XYZ;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }

	    if(bc->bcX1 == BC::SPONGE){
		FOR_X{
		    if(domain->x[i] > domain->Lx - spongeLX){
		        double spongeX = (domain->x[i] - (domain->Lx - spongeLX))/spongeLX;
		        FOR_Y{
			    FOR_Z{
			        int ii = GET3DINDEX_XYZ;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeX, 2.0) + 0.845*pow(spongeX, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }    

	    if(bc->bcY0 == BC::SPONGE){
		FOR_Y{
		    if(domain->y[j] < spongeLY){
		        double spongeY = (spongeLY - domain->y[j])/spongeLY;
		        FOR_X{
			    FOR_Z{
			        int ii = GET3DINDEX_XYZ;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }
	
	    if(bc->bcY1 == BC::SPONGE){
		FOR_Y{
		    if(domain->y[j] > domain->Ly - spongeLY){
		        double spongeY = (domain->y[j] - (domain->Ly - spongeLY))/spongeLY;
		        FOR_X{
			    FOR_Z{
			        int ii = GET3DINDEX_XYZ;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeY, 2.0) + 0.845*pow(spongeY, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }    

	    if(bc->bcZ0 == BC::SPONGE){
		FOR_Z{
		    if(domain->z[k] < spongeLZ){
		        double spongeZ = (spongeLZ - domain->z[k])/spongeLZ;
		        FOR_X{
			    FOR_Y{
			        int ii = GET3DINDEX_XYZ;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }
	
	    if(bc->bcZ1 == BC::SPONGE){
		FOR_Z{
		    if(domain->z[k] > domain->Lz - spongeLZ){
		        double spongeZ = (domain->z[k] - (domain->Lz - spongeLZ))/spongeLZ;
		        FOR_X{
			    FOR_Y{
			        int ii = GET3DINDEX_XYZ;
			        sigma[ii] = fmax(spongeStrength*(0.068*pow(spongeZ, 2.0) + 0.845*pow(spongeZ, 8.0)), sigma[ii]);
			    }
		        }
		    }
		}
	    }    



	    std::cout << " > Done initializing sponge!" << std::endl;

	}

};

#endif
