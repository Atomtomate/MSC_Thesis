#ifndef SOHUBBARD_H_
#define SOHUBBARD_H_

#include "DMFT_BetheLattice.hpp"
#include "DMFTSolver.hpp"
#include "Config.hpp"
#include "integrals.hpp"

#include "Faddeeva.hh"
#include "gnuplot-iostream.h"

#include <iostream>
#include <string>
#include <thread>

namespace DMFT{
namespace examples{


   int _test_SOH(DMFT::Config& config, bool use_bethe, double mxing);
   void _run_hysteresis(double);
   void _test_hysteresis(void);


}	//end namespace examples
}   //end namespace DMFT
#endif
