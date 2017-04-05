#ifndef CONFIG_H_
#define CONFIG_H_
// ========== set up logging ==========
#define ELPP_STL_LOGGING
//#define ELPP_DISABLE_VERBOSE_LOGS
//#define ELPP_DISABLE_DEBUG_LOGS
#include "easylogging++.h"
//ELPP_DISABLE_LOGS

#include "gnuplot-iostream.h"

#include "mpi.h"
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

//#define TIME_ORDERED_CONFIG_LIST

#include <utility>		// pair
#include <algorithm>    	// std::binary_search, std::sort
#include <set>			// multiset
#include <vector>
#include <cstddef>		// nullptr

#include <boost/cstdfloat.hpp>
#include <boost/math/constants/constants.hpp>

//#define EIGEN_NO_DEBUG // only for release build
//#define EIGEN_USE_MKL_ALL //huge overhead for user time (1 cpu -> 2 cpus)

#include <Eigen/Dense>


#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>


namespace DMFT
{

    // ========== global config parameters ==========
    // TODO: get this from user, define static config to be passed around
    // TODO: consider adjustable maxMatsFreq
    // TODO: export some of this for weakCoupling only
    static const	int _CONFIG_maxMatsFreq = 2048;//4096;		// frequencies go from maxMatsFreq/2 to maxMatsFreq/2 - 1
    static const    int _CONTIG_minMF    = -static_cast<int>(_CONFIG_maxMatsFreq/2.0);
    static const    int _CONFIG_maxTBins =	2048;//8192;//4096; 		// (default default 131072, 65536) powers of 2 make fft faster, both bins need to be equal for FFT
    static const	int _CONFIG_maxSBins = 16384;			// TODO: dependent on tBins, this is acutalle sBinsSize/beta
    // TODO: assert that spins are integeres starting with 0 (indexing!)
    static const    int _CONFIG_spins = 2;
    enum SPIN {DOWN = 0, UP = 1};

    //#define DEBUG_MODE


    // ========== Typedefs for abstraction from underlying libraries ==========
    typedef double RealT;			// TODO: better adjustment for precision
    typedef std::complex<RealT> ComplexT;		// this is PETSC complex, when compiled with --with-clanguage=c++


    // ========== Matrix typedefs. For now this aliases Eigen3 ==========
    typedef Eigen::MatrixXd	    MatrixT;
    typedef Eigen::RowVectorXd  RowVectorT;
    typedef Eigen::VectorXd	    VectorT;
    typedef Eigen::MatrixXcd    CMatrixT;
    typedef Eigen::VectorXcd    CVectorT;
    typedef Eigen::RowVectorXcd CRowVectorT;

    //TODO: make length dynamic for larger malues of _CONFIG_maxMatsFreq and use config constructor to intiialize
    typedef Eigen::ArrayXXcd MatG;
    //TODO: assert maxMatsFreq == maxTBins
    typedef Eigen::ArrayXXd ImTG;


    // ========== Additional typedefs ==========
    typedef Eigen::ArrayXXd Potential;
    typedef Eigen::ArrayXd Energies;
    typedef Eigen::ArrayXcd Hybridization;

    static const RealT PI = boost::math::constants::pi<RealT>();

    typedef Eigen::Matrix<RealT, _CONFIG_maxMatsFreq, 1> FixedMFGrid;
    typedef VectorT MFGrid;
    // 2 pi (n+1)/ beta
    //typedef Mat MatrixT;


    // ========== Utility Functions ==========

    template <typename T>
        inline int sgn(const T val) {return (T(0)< val) - (val < T(0));}

    inline RealT mRInd(const int freq, const RealT beta) {return (2*freq+1)*PI/beta;}

    //#define mFreq(freq,beta) (2.0*((freq)-static_cast<int>(_CONFIG_maxMatsFreq/2.0))+1)*boost::math::constants::pi<RealT>()/(beta)
    inline RealT mFreq(const int freq,const RealT beta) {
        return (2.0*(freq-static_cast<int>(_CONFIG_maxMatsFreq/2.0))+1)*boost::math::constants::pi<RealT>()/beta;
    }

    inline long int arrIndex(const int i, const int spin) {return spin+i*_CONFIG_spins;}


    // ========== Config Struct ==========
    /*! This stores all high level parameters for the DMFT solver.
     *  
     */
    struct Config
    {
        public:
            // TODO: decide whether to read config from file or use macros
            /*! 
             *  param  [in]  beta       inverse temperature
             *  param  [in]  mu         chemical potential
             *  param  [in]  U          interaction strength
             *  param  [in]  mfCount    number of Matsubara frequencies, a grid with all frequencies will be cosntructed, max floor((mfCount-1)/2) 
             */
            Config(const RealT beta,const RealT mu,const RealT U,const int mfCount): beta(beta), mu(mu), mfCount(mfCount), mfGrid(mfCount) ,U(U)
        {
            const int min = static_cast<int>(mfCount/2.0);
            for(int n= -static_cast<int>(mfCount/2.0); n<static_cast<int>((mfCount-1)/2.0);n+=1)
            {
                mfGrid(n + min) = PI*(2.0*n+1)/beta;
            }
        }

            //TODO: sparse mfGrid with interpolation

            const RealT beta;
            const RealT mu;
            const RealT U;
            const int mfCount;
            MFGrid mfGrid;

    };


}   //end namespace DMFT
#endif
