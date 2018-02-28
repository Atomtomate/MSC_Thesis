#ifndef IMPSOLV_HPP_
#define IMPSOLV_HPP_
#include "Config.hpp"
#include "GreensFct.hpp"

namespace DMFT
{



/*!	This is the base class for all impurity Solvers.
 *	Using CRTP (static polymorphism) we can avoid calls to virtual functions.
 */
template <class Implementation>
class ImpSolver
{
	public:
		inline void update(const unsigned long iterations = 10000l)     const { 	this->impl().update(iterations); }
		inline void computeImpGF(void) 					const { 	this->impl().computeImpGF(); }
		inline GreensFct *const getImpGF(void) 				const { return	this->impl().getImpGF(); }
        inline GreensFct *const getWeissGF(void)                              const { return	this->impl().getWeissGF(); }
        inline int expansionOrder(void)                                 const { return  this->impl().expansionOrder(); }
        inline RealT avgN(void)                                         const { return  this->impl().avgN(); }

    protected:

	private:
		inline Implementation& impl()
		{
			return *static_cast<Implementation*>(this);
		}
};

}
#endif
