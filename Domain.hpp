#ifndef _DOMAINH_
#define _DOMAINH_

#include <iostream>

class Domain{

    public:

        int Nx, Ny, Nz, N;
        double *x, *y, *z;
        double Lx, Ly, Lz;
	double dx, dy, dz;


    Domain(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz){
	
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	N = Nx*Ny*Nz;

	this->Lx = Lx;
	this->Ly = Ly;
	this->Lz = Lz;

	x = new double[Nx];	
	y = new double[Ny];	
	z = new double[Nz];

        for(int ip = 0; ip < Nx; ip++){
	    x[ip] = (((double)ip)/((double)Nx - 1.0))*Lx;
	}	

        for(int jp = 0; jp < Ny; jp++){
	    y[jp] = (((double)jp)/((double)Ny - 1.0))*Ly;
	}	

        for(int kp = 0; kp < Nz; kp++){
	    z[kp] = (((double)kp)/((double)Nz - 1.0))*Lz;
	}	

	dx = x[1]-x[0];
	dy = y[1]-y[0];
	dz = z[1]-z[0];

	std::cout << " >Domain initialization..." << std::endl;
	std::cout << " >Domain: " << Lx << "x" << Ly << "x" << Lz << std::endl;
	std::cout << " >Mesh: " << Nx << "x" << Ny << "x" << Nz << std::endl;
	std::cout << " >Total Points: " << Nx*Ny*Nz << std::endl; 

    }

};

#endif
