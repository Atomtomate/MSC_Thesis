#ifndef MAIN_H_
#define MAIN_H_

//#include <mgl2/mgl.h>
//#include <mgl2/qt.h>

#include "Config.hpp"
#include "singleOrbHubbard.hpp"

#include <iostream>
#include <fstream>
#include <string>


#include <boost/mpi.hpp>

//#include <petscsys.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

// numeric integration
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::accumulators;
using boost::accumulators::accumulator_set;
using boost::accumulators::features;

//#include <boost/mpi/environment.hpp>
//#include <boost/mpi/communicator.hpp>
//namespace mpi = boost::mpi;

// boost tends to throw compilation bugs when placing IOHelper include before boost include


//void _TEST_fftw3(void);

#endif
