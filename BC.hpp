#ifndef _BCH_
#define _BCH_

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

	}

};

#endif
