#ifndef LATTICE_HPP_ #define LATTICE_HPP_

#include "Config.hpp"

namespace DMFT
{


/*!	This is the base class for all lattice types.
 *	Using CRTP (static polymorphism) we can avoid calls to virtual functions.
 */
template <class Implementation>
class Lattice
{
	public:
	private:
		inline Implementation& impl()
		{
			return *static_cast<Implementation*>(this);
		}
};

/*! Bethe lattice with half-bandwidth D.
 */
class BetheLattice: public Lattice<BetheLattice>
{
    public:
        BetheLattice(const RealT D): D(D) {}

        /*! Compute and return D(e) for a semi circular DOS.
         *  \f[ D( \epsilon ) = \frac{ \sqrt{4 t^2 - \epsilon^2} }{2 \pi t^2 } | \epsilon | < 2t  \f] 
         *  TODO: implement
         */
        inline RealT semiCircDOS(RealT eps)
        {
            return 0;
        }

        const RealT D;  // half bandwidth
};

}   //end namepsace DMFT

#endif
