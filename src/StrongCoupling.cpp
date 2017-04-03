#include "StrongCoupling.h"
#include <iostream>


namespace  DMFT
{
/*StrongCoupling::StrongCoupling(SConfigL& confs, MatrixT& hybr, const Realt U, const RealT zeroShift, const RealT beta, const unsigned int burninSteps):
        hybr(hybr), confs(confs), U(U), zeroShift(zeroShift), beta(beta), burninSteps(burninSteps), g0(I, beta)
{
        auto seed = 1234567890;
        // auto seed = make_sprng_seed();
        init_sprng(seed,SPRNG_DEFAULT, DEFAULT_RNG_TYPE);
        this->n = hybr.rows();
        steps = 0;
        //g0(maxMatsFreq,spins);          //TODO: impl. hybr
}*/

StrongCoupling::~StrongCoupling()
{
        //ptr->free_sprng();                                                    //only needed for default interface
}

void StrongCoupling::doMeasurement(void)
{
	
	//(M.transpose() * M2).diagonal(); // m[i,:].dot(m[:,i])
}

RealT StrongCoupling::acceptanceR(const RealT U,const RealT beta) const
{
        // TODO: placeholder
        return 0.0;
}

int StrongCoupling::update(const RealT U,const RealT beta)
{
	auto	zeta	= 0;
        auto    zetap   = 1;
        char    sigma   = 1;
        unsigned int n  = hybr.rows();
        steps           += 1;
	
	if(zeta < 0.5){							// try to insert/remove segment
		if(zetap < 0.5){					// try to insert segment
		
		} else {						// try to remove segment

		}
	} else {							// try to insert/remove anti-segment
		if(zetap < 0.5){					// try to insert anti-segment
		
		} else {						// try to remove anti-segment

		}
	}
return 0;	
}

} //end namespace DMFT
