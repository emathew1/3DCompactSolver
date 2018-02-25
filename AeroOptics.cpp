#include "AeroOptics.hpp"

void AeroOptics::generateBaseBeamDiff(int angleNumber){


    //We know our y0 and y1 from our upper and lower IBound's
    double y0  = lowerIBound;
    double y1  = upperIBound;

    double minAngle = 0.0;
    double currentAngle = minAngle + ((double)angleNumber/((double)(numAngles-1)))*(maxAngle-minAngle);
    currentAngle *= M_PI/180.0;

    //Get the length of the integration line for this angle
    Ldiag = totalDY/cos(currentAngle);
    totalDX = Ldiag*sin(currentAngle);

    //Number of points for this line
    Nline  = ceil(Ldiag/ds) + 1;
    dsLine = Ldiag/(double)(Nline - 1);
    dxLine = totalDX/(double)(Nline - 1);
    dyLine = totalDY/(double)(Nline - 1); 
}

void AeroOptics::writeOPLFiles(){


    //Going thru all of the angles
    for(int ap = 0; ap < numAngles; ap++){
  
	double minAngle = 0.0;  
	double currentAngle = minAngle + ((double)ap/((double)(numAngles-1)))*(maxAngle-minAngle);
	char buff[64];
	sprintf(buff, "%4.2f", currentAngle); 
	string filename = buff;
	filename = filename + ".OPL";

	ofstream fp(filename, ios::binary | ios::app);
	fp.write((char*)&OPL[ap*Nx*Nz], sizeof(double)*Nx*Nz);   
	fp.close();
    }

}

void AeroOptics::computeAO(){

    //Going thru all of the angles...
    for(int ap = 0; ap < numAngles; ap++){

	//generate the geometry we need to construct the ray line...
	generateBaseBeamDiff(ap);

	//Going thru all of the positions in the x-z plane
        FOR_Z{
	    FOR_X{
		
		//this is our integration "ray"
	        double *beamRho = new double[Nline];
			
		//interpolating the density along the line...
	        for(int sp = 0; sp < Nline; sp++){
		    double xp_temp[3];

		    //We can come out the other side in a periodic domain in x 
		    xp_temp[0] = fmod(domain->x[i]  + (double)(sp*dxLine),domain->Lx);
		    xp_temp[1] = lowerIBound   + (double)(sp*dyLine);
		    xp_temp[2] = domain->z[k]; 

		    beamRho[sp] = linearInterpolation(domain, rho, xp_temp);
	        }

		//take the integral along the beam direction...
		//trapezoidal rule...
		double localOPL = beamRho[0]*dsLine/2.0 + beamRho[Nline-1]*dsLine/2.0;
		for(int sp = 1; sp < Nline-1; sp++)
		    localOPL += beamRho[sp]*dsLine;

		OPL[ap*Nx*Nz+k*Nx+i] = localOPL; 
					//This is really only a part of the OPL
					//OPL = \int_0^L n ds = \int_0^L (1 + K_GD \rho) \ds
					//OPL = L + K_GD *  [\int_0^L rho ds]
					//		    ^
					//		    |--------This part
	        delete[] beamRho;		   
	    }
         }
    }

    //write a function to output the OPL files...

}


