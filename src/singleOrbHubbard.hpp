#ifndef SOHUBBARD_H_
#define SOHUBBARD_H_

#include "DMFT_BetheLattice.hpp"
#include "DMFTSolver.hpp"
#include "Config.hpp"
#include "integrals.hpp"
#include "FFT.hpp"

#include "Faddeeva.hh"
#include "gnuplot-iostream.h"

#include <iostream>
#include <string>
#include <thread>

namespace DMFT{
namespace examples{


   int _test_PT(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator, bool use_bethe, double mxing);
   int _test_SOH(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator);
   void _run_hysteresis(RealT beta, boost::mpi::communicator world, boost::mpi::communicator local, const bool isGenerator);
   void _test_hysteresis(boost::mpi::communicator world, boost::mpi::communicator local, const bool isGenerator);


}	//end namespace examples
}   //end namespace DMFT
#endif
