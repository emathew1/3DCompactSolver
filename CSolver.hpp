#include "Macros.hpp"
#include "Utils.hpp"
#include "BC.hpp"
#include "TimeStepping.hpp"
#include "IdealGas.hpp"
#include "SpongeBC.hpp"

class CSolver{

    public:

	Domain *domain;
	BC *bc;
	TimeStepping *timeStepping;
	IdealGas *idealGas;
	double alphaF;
	double mu_ref;

	CSolver(Domain *domain, BC *bc, TimeStepping *timeStepping, double alphaF, double mu_ref){

	    this->domain = domain;
	    this->bc = bc;
	    this->timeStepping = timeStepping;
	    this->alphaF = alphaF;
	    this->mu_ref = mu_ref;

	    idealGas = new IdealGas(domain);

	}

};
