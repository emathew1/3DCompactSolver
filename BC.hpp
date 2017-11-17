#ifndef _BCH_
#define _BCH_

#include <iostream>

class BC{

    public:

	enum BCType {PERIODIC_SOLVE, DIRICHLET_SOLVE};
	enum BCKind {PERIODIC, SPONGE, ADIABATIC_WALL, MOVING_WALL, INLET};

	BCType bcXType, bcYType, bcZType;
	BCKind bcX0, bcX1, bcY0, bcY1, bcZ0, bcZ1;

	BC(BCType bcXType, BCKind bcX0, BCKind bcX1, 
	   BCType bcYType, BCKind bcY0, BCKind bcY1,
	   BCType bcZType, BCKind bcZ0, BCKind bcZ1){

	    this->bcXType = bcXType;
	    this->bcYType = bcYType;
	    this->bcZType = bcZType;

	    this->bcX0 = bcX0;
	    this->bcX1 = bcX1;
	    this->bcY0 = bcY0;
	    this->bcY1 = bcY1;
	    this->bcZ0 = bcZ0;
	    this->bcZ1 = bcZ1;

	    std::cout << std::endl;
	    std::cout << " > Initializing boundary conditions..." << std::endl;
	
	    std::cout << " >     X BOUNDARY CONDITIONS    " << std::endl;
	    if(bcXType == PERIODIC_SOLVE){
		std::cout << " > ----------PERIODIC---------- " << std::endl;
	    }else{
		if(bcX0 == SPONGE){
		    std::cout << " > X0=SPG";
		}else if(bcX0 == ADIABATIC_WALL){
		    std::cout << " > X0=ABW";
		}else if(bcX0 == MOVING_WALL){
		    std::cout << " > X0=MOW";
		}

		if(bcX1 == SPONGE){
		    std::cout << "----------------SPG=X1" << std::endl;
		}else if(bcX1 == ADIABATIC_WALL){
		    std::cout << "----------------ABW=X1" << std::endl;
		}else if(bcX1 == MOVING_WALL){
		    std::cout << "----------------MOW=X1" << std::endl;
		}

	    }

	    std::cout << " >     Y BOUNDARY CONDITIONS    " << std::endl;
	    if(bcYType == PERIODIC_SOLVE){
		std::cout << " > ----------PERIODIC---------- " << std::endl;
	    }else{
		if(bcY0 == SPONGE){
		    std::cout << " > Y0=SPG";
		}else if(bcY0 == ADIABATIC_WALL){
		    std::cout << " > Y0=ABW";
		}else if(bcY0 == MOVING_WALL){
		    std::cout << " > Y0=MOW";
		}

		if(bcY1 == SPONGE){
		    std::cout << "----------------SPG=Y1" << std::endl;
		}else if(bcY1 == ADIABATIC_WALL){
		    std::cout << "----------------ABW=Y1" << std::endl;
		}else if(bcY1 == MOVING_WALL){
		    std::cout << "----------------MOW=Y1" << std::endl;
		}

	    }

	    std::cout << " >     Z BOUNDARY CONDITIONS    " << std::endl;
	    if(bcZType == PERIODIC_SOLVE){
		std::cout << " > ----------PERIODIC---------- " << std::endl;
	    }else{
		if(bcZ0 == SPONGE){
		    std::cout << " > Z0=SPG";
		}else if(bcZ0 == ADIABATIC_WALL){
		    std::cout << " > Z0=ABW";
		}else if(bcZ0 == MOVING_WALL){
		    std::cout << " > Z0=MOW";
		}

		if(bcZ1 == SPONGE){
		    std::cout << "----------------SPG=Z1" << std::endl;
		}else if(bcZ1 == ADIABATIC_WALL){
		    std::cout << "----------------ABW=Z1" << std::endl;
		}else if(bcZ1 == MOVING_WALL){
		    std::cout << "----------------MOW=Z1" << std::endl;
		}
	     }
	

	}

};

#endif
