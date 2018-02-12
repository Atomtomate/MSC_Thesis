#ifndef STRONG_COUPLING_HPP_
#define STRONG_COUPLING_HPP_


#include "Config.hpp"
#include "GreensFct.hpp"
#include "ImpSolver.hpp"
#include "IOhelper.hpp"
#include "GFLPoly.hpp"
#include "Segments.hpp"

#include <boost/serialization/vector.hpp>
#include <boost/accumulators/numeric/functional.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <iostream>
#include <tuple>
#include <algorithm>
#include <limits>
//#include <petscsys.h>                         // random numbers

#define _CT_HYB_USE_L_BASIS 0
namespace DMFT
{
/** Impurity Sovler using the CT-HYB algorithm.
 *  The following data structures are used: \n
 *  Configurations:     \f$C = \{(\tau_1_up,taup_1_up),(\tau_2_up,taup_2_up),...(\tau_1_down,taup_1_down),(\tau_2_down,taup_2_down),...\}\f$ \n
 *  Hybridization function: \f$Hybr\f$ \n
 *  Inverse Hybridization function matrix M:
 *      \f$[M^{-1}_{\sigma}]_{ij} = G^{\sigma}_0(\tau_i - \tau_j) - \alpha_{\sigma}(s_i) \delta_{ij}\f$ \n
 *  Interaction: U \n
 *  Inverse Temperature: \f$\beta\f$ \n
*/



class StrongCoupling;
typedef void (StrongCoupling::*MC_Move)(void);
struct handler_pair
{
    int code;
    MC_Move fn;
};

class StrongCoupling
{
    public:
        using AccT = boost::accumulators::accumulator_set<RealT, boost::accumulators::features<boost::accumulators::tag::sum, boost::accumulators::tag::variance > >;
        using AccMT = boost::accumulators::accumulator_set<RealT, boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::skewness, boost::accumulators::tag::kurtosis > >;
        using AccCT = boost::accumulators::accumulator_set<ComplexT, boost::accumulators::features<boost::accumulators::tag::sum, boost::accumulators::tag::moment<2> > >;
        StrongCoupling(GreensFct* const hybr, GreensFct* const gImp, const Config& config, const unsigned int burninSteps);
        virtual ~StrongCoupling();


        /**     Updates the time ordered spin configuration and the inverse
         *  weiss greens function M by inserting or removing one configuration.
         *  @param C List of configurations
         *  @param M Inverse Weiss Green's function
         *  @param G0 Weiss Green's function
         *  @param U Interaction
         *  @param beta Inverse temperature
         */
        int update(const unsigned long int iterations = 1l);
        RealT avgN(void)
        {
            RealT res = 0;
            for(int f =0; f < _CONFIG_spins; f++)
            {
                res += boost::accumulators::mean(expOrd[f]);
            }
            return res/_CONFIG_spins;
        }
        int expansionOrder(void)
        {
            if(steps>burninSteps)
                return boost::accumulators::mean(expOrd[0]);
            return 0;
        }

        inline GreensFct * const getImpGF(void) {
            return gImp;
        }
        void computeImpGF(void);


        /** computes the acceptance rate by evaluation of the determinant
         * ratios \f$ \frac{-\beta U}{(n+1}  \prod \limits_{\sigma}
         * \frac{det[M(n+1,\sigma)^{-1}]}{det[M(n,\sigma)^{-1}]}\f \quad GKW (8.36)$
         */
        RealT acceptanceR(const RealT U, const RealT beta) const;

    private:
        using ProposalRes = std::pair<RealT, bool>;
        trng::yarn2 r_time, r_timep, r_spin, r_insert, r_accept, r_shift;   // random number engines
        trng::uniform01_dist<> u;                                           // random number distribution
        GreensFct* const hyb;
        GreensFct* const gImp;
        GFLPoly gImpLPoly;
        const Config &conf;
        std::array<MatrixT,_CONFIG_spins> M;
        Segments<_CONFIG_spins> segments;
        // segmentCache[FLAVOR][TIME_OF_INSERT]. {first - segment start time, second - segment end time}
        std::array<std::vector<std::pair<RealT,RealT>>, _CONFIG_spins> segmentCache;

        unsigned int steps;                                             // number of updates
        const unsigned int burninSteps;

        int lastSign;						                            // needed when proposal is rejected
        long int totalSign;
        std::array<AccMT, _CONFIG_spins> expOrd;
        std::array<std::array< AccT, _CONFIG_maxTBins>, _CONFIG_spins> itBins;
        std::array<std::array<ComplexT, _CONFIG_maxMatsFreq>, _CONFIG_spins> mfBins;
        std::array<std::array<RealT,_CONFIG_maxLPoly>,_CONFIG_spins> gl_c;

        decltype(auto) tryZeroToFull(const unsigned f, const RealT fac, const RealT zetap);
        decltype(auto) tryFullToZero(const unsigned f, const RealT fac, const RealT zetap);
        decltype(auto) tryInc(const RealT t, const RealT tp, const unsigned int f, const RealT fac, const RealT zetap = -1.);
        decltype(auto) tryDec(const unsigned int row, const unsigned int f, const RealT fac, const RealT zetap = -1.);
        decltype(auto) tryInsAntiSeg(const RealT t_n, const RealT tp_n, const unsigned int f_n, const RealT fac, const RealT zetap = -1.);
        decltype(auto) tryRemAntiSeg(int index, const unsigned int f_n, const RealT fac, const RealT zetap = -1.);

        void updateContribution(void);
        void updateLPoly(void);


        /*! Positions row and col at index from at index to
         *  @param [out] A Matrix to be transformed
         *  @param [in] from row/col index
         *  @param [in] to   target row/col index
         */
        void swapRows(MatrixT *A, const int from, const int to);
        void MInc(MatrixT *A, const RowVectorT &R, const VectorT &Q, const RealT Sp, const int index);
        void MDec(MatrixT* A, const unsigned int row);

        inline RealT hybCall(RealT ts, RealT tf, const unsigned char flavor)
        {
            RealT t = std::fmod(tf,conf.beta) - std::fmod(ts,conf.beta);
            //return -hyb.getByT(t, flavor, 1);
            return -(*hyb)(t,flavor);
        }

        // for debugging purposes
        std::array<RealT,_CONFIG_spins> fullDet;
        std::array<RealT,_CONFIG_spins> runningDet;
        std::array<MatrixT,_CONFIG_spins> Minv;
        void debug_test(const unsigned int f_n);
        void rebuildM(const bool timeOrdered = false);
};

} //end namespace DMFT
#endif

