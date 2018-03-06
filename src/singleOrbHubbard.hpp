#ifndef SOHUBBARD_H_
#define SOHUBBARD_H_

#include "DMFT_BetheLattice.hpp"
#include "IPT.hpp"
#include "DMFTSolver.hpp"
#include "Config.hpp"
#include "integrals.hpp"
#include "FFT.hpp"
#include "GFTail.hpp"
#include "Segments.hpp"

#include "Faddeeva.hh"
#include "gnuplot-iostream.h"

#include <iostream>
#include <string>
#include <thread>

namespace DMFT{
namespace examples{


   int _test_hyb(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator);
   void _test_IPT(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator);
   void _IPT_PD(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator);
   void _IPT_Z(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator);
   int _test_PT(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator, bool use_bethe, double mxing);
   int _test_SOH(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator);
   void _test_average_PO(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator);
   void _run_hysteresis(RealT beta, boost::mpi::communicator world, boost::mpi::communicator local, const bool isGenerator);
   void _test_hysteresis(boost::mpi::communicator world, boost::mpi::communicator local, const bool isGenerator);


}	//end namespace examples
}   //end namespace DMFT
#endif
