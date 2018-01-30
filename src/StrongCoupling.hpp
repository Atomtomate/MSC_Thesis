#ifndef STRONG_COUPLING_HPP_
#define STRONG_COUPLING_HPP_


#include "Config.hpp"
#include "GreensFct.hpp"
#include "ImpSolver.hpp"
#include "Segments.hpp"

#include <boost/serialization/vector.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

#include <iostream>
#include <tuple>
#include <algorithm>
#include <limits>
//#include <petscsys.h>                         // random numbers


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





class StrongCoupling
{
    public:
        using AccT = boost::accumulators::accumulator_set<RealT, boost::accumulators::features<boost::accumulators::tag::sum, boost::accumulators::tag::moment<2> > >;
        StrongCoupling(GreensFct &hybr, GreensFct &gImp, const Config& config, const unsigned int burninSteps);
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
        RealT avgN(void) {return -1;}
        int expansionOrder(void) {return -1;}
        void computeImpGF(void) {};
        void doMeasurement(void);


        /** computes the acceptance rate by evaluation of the determinant
         * ratios \f$ \frac{-\beta U}{(n+1}  \prod \limits_{\sigma}
         * \frac{det[M(n+1,\sigma)^{-1}]}{det[M(n,\sigma)^{-1}]}\f \quad GKW (8.36)$
         */
        RealT acceptanceR(const RealT U, const RealT beta) const;

    private:
        trng::yarn2 r_time, r_timep, r_spin, r_insert, r_accept, r_shift;   // random number engines
        trng::uniform01_dist<> u;                                           // random number distribution
        GreensFct &hyb;
        GreensFct &gImp;
        const Config &conf;
        std::array<MatrixT,_CONFIG_spins> M;
        Segments<_CONFIG_spins> segments;
        // segmentCache[FLAVOR][TIME_OF_INSERT]. {first - segment start time, second - segment end time}
        std::array<std::vector<std::pair<RealT,RealT>>, _CONFIG_spins> segmentCache;

        unsigned int steps;                                             // number of updates
        const unsigned int burninSteps;

        int lastSign;						                            // needed when proposal is rejected
        long totalSign;
        long totN;
        std::array< AccT, _CONFIG_maxSBins> itBinsUP;
        std::array< AccT, _CONFIG_maxSBins> itBinsDOWN;

        RealT tryInc(const RealT t, const RealT tp, const unsigned char f, const RealT fac, const RealT zetap = -1.);
        RealT tryDec(const unsigned int row, const unsigned char f, const RealT fac, const RealT zetap = -1.);
        RealT MInc(MatrixT *A, const RowVectorT &R, const VectorT &Q, const RealT Sp, const int index);
        void MDec(MatrixT* A, const unsigned int row);
        bool tryInsAntiSeg(const RealT t_n, const RealT tp_n, const int f_n, const RealT fac, const RealT zetap = -1.);
        bool tryRemAntiSeg(int index, const int f_n, const RealT fac, const RealT zetap = -1.);

        inline RealT hybCall(RealT ts, RealT tf, const unsigned char flavor)
        {
            if(tf > conf.beta) tf -= conf.beta;
            const RealT t = tf - ts;
            return hyb(t,flavor);
        }
};

} //end namespace DMFT
#endif

