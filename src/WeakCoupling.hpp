#ifndef WEAK_COUPLING_H_
#define WEAK_COUPLING_H_

#include "Config.hpp"
#include "GreensFct.hpp"
#include "ImpSolver.hpp"
#include "MCAccumulator.hpp"

#include <boost/serialization/vector.hpp>

#include <iostream>
#include <tuple>
#include <algorithm>

#define MEASUREMENT_SHIFT 1
#define MATSUBARA_MEASUREMENT 1

namespace DMFT
{

    //#define WK_NAIVE_MEASUREMENT

    /*! collection of functions for the impurity solver using the CT-INT algorithm.
     *  The following data structures are used: \n
     *  Configurations:	\f$C = \{(\tau_1,s_1),(\tau_2,s_2),...\}\f$ \n
     *  Weiss green's function: \f$G_0\f$ \n
     *  Inverse Weiss Green's function matrix M:
     *  	\f$[M^{-1}_{\sigma}]_{ij} = G^{\sigma}_0(\tau_i - \tau_j) - \alpha_{\sigma}(s_i) \delta_{ij}\f$ \n
     *  Interaction: U \n
     * 	Inverse Temperature: \f$\beta\f$ \n
     */

    class WeakCoupling: public ImpSolver<WeakCoupling>
    {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

                // ========== Definitions ==========
                /*! Constructor for the WeakCoupling solver
                 *  @param [in]  g0				Weiss Green's function
                 *  @param [out] gImp			sampled impurity Green's function
                 *  @param [in]	 U				Interaction strength
                 *  @param [in]  zeroShift		U-zeroShift avoids some of the sign problem complications
                 *  @param [in]  mu				chemical potential
                 *  @param [in]  beta			temperature
                 *  @param [in]  burningSteps	disregard first values during simulation until reasonable stable state is achieved TODO: auto gen
                 */
                WeakCoupling(GreensFct &g0, GreensFct &gImp, const Config& config, const RealT zeroShift, const unsigned int burninSteps);
            virtual ~WeakCoupling();

            inline int expansionOrder(void) {return n;}

            /*!	@brief	Updates the time ordered spin configuration and the inverse
             *  		Weiss Greens function M by inserting or removing one configuration.
             *
             *	@param [in]	iterations	number of updates. If left empty: one update
             */
            void update(const unsigned long iterations = 1l);

            /*!	@brief	computes imaginary time impurity Green's function after sampling
             *
             *	After each update, the helper function updateContribution is called, which
             *	collects the contribution of the current confiuration into bins, these are
             *	measured later. See Gull et al. (197)
             *
             *	@return	measured iTime Green's function 
             */
            void computeImTGF(void);

            /*!	@brief	computes Matsubara impurity Green's function after sampling
             *		
             *	This uses the same binned values obtained during sampling (see: updateContribuiton)
             *	as computeImTGreensFct. Since complex exponentials are needed, there is some additional
             *	computational effort. See Gull et al. (198) 
             * 	Better results for high frequencies, if not needed, meassure iTime and FFT
             * 	
             *	@return	measured Mats Green's function 
             */
            void computeMatGreensFct(void);

            /*!	@brief	comppute both imaginary time and Matsubara impurity GF.
             *	@return	impurity GF
             */
            void computeImpGF(void);


            void computeImpGF_OLD(void);
            /*! @return impurity Green's function. sample and call compute first!
            */
            inline GreensFct& getImpGF(void) {

                //if(gImp_NeedsUpdate) LOG(WARNING) << "Impurity Green's function requested but not compted";
                return gImp;
            }

            inline GreensFct& getWeissGF(void) {

                //if(gImp_NeedsUpdate) LOG(WARNING) << "Impurity Green's function requested but not compted";
                return g0;
            }

            /* ratios \f$ \frac{-\beta U}{(n+1}  \prod \limits_{\sigma}
             * \frac{det[M(n+1,\sigma)^{-1}]}{det[M(n,\sigma)^{-1}]}\f \quad GKW (8.36)$
             */
            void acceptanceR(void);

            /*! rebuild M_\sigma from scratch by direct inversion
             *
             *  @param [in] spin \sigma
             *
             *  @result M_\sigma
             */
            MatrixT rebuildM(int spin);

        private:
            // ========== Typedef for <spin,itime> configurations ==========
            typedef std::tuple<double, int> SConfig;						// time, spin, sign
            struct compare {
                inline bool operator()(const SConfig& lhs, const SConfig& rhs) const {
                    return std::get<0>(lhs) < std::get<0>(rhs);}
            };
            typedef std::set<SConfig, compare> SConfigTOL;  				// Time orderer list (set)
            typedef std::vector<SConfig> SConfigSOL; 						// Ordered by time of generation
#ifdef TIME_ORDERED_CONFIG_LIST
            typedef SConfigTOL SConfigL; 								// TODO: switch to Boost::set
            inline void pushConfig(const SConfig& c){
                confs.insert(c);
            }
            inline void deleteConfig(const int pos){
                auto sConfIt = confs.begin();
                std::advance(sConfIt, pos);
            }
#else															// TODO: reseve size
            typedef SConfigSOL SConfigL;
            inline void pushConfig(const SConfig& c){
                confs.push_back(c);
                const auto gfSize = gfCache[UP].size();
                gfCache[UP].conservativeResize(gfSize+1);
                gfCache[UP](gfSize) = g0(std::get<0>(c),UP);
                gfCache[DOWN].conservativeResize(gfSize+1);
                gfCache[DOWN](gfSize) = g0(std::get<0>(c),DOWN);
            }
            //TODO: erase requires non const lvalue, implement cache
            inline void deleteConfig(const int pos){
                confs.erase(confs.begin() + pos);
                LOG(WARNING) << "delete cache not implemented yet";
            }
            inline void popConfig(void){
                confs.pop_back();
                const auto gfSize = gfCache[UP].size();
                gfCache[UP].conservativeResize(gfSize-1);
                gfCache[DOWN].conservativeResize(gfSize-1);
            }

            inline void swapConfigs(const int pos1, const int pos2){
                std::swap(confs[pos1], confs[pos2]);
                RealT tmp;
                tmp = gfCache[UP](pos1);
                gfCache[UP](pos1) = gfCache[UP](pos2);
                gfCache[UP](pos2) = tmp;
                tmp = gfCache[DOWN](pos1);
                gfCache[DOWN](pos1) = gfCache[DOWN](pos2);
                gfCache[DOWN](pos2) = tmp;
            }
#endif
            trng::yarn2 r_time, r_spin, r_insert, r_accept, r_shift; // random number engines
            trng::uniform01_dist<> u; // random number distribution

            //REMARK: there should probably be some check of ownership, smart pointer could be used but are slow
            GreensFct   &g0;			// reference to Weiss GF
            GreensFct   &gImp;			// sampled impurity GF

            const Config& config;
            std::array<MatrixT,2> M;
            std::array<VectorT,2> gfCache;	                // cached access to Weiss GF at current vertex points
            SConfigL confs;

            unsigned int steps;					// number of updates
            const unsigned int burninSteps;		// throw away some steps at the start
            int lastSign;						// needed when proposal is rejected
            int totalSign;
            int n; 					// expansion order (number of used rows/cols)
            const RealT zeroShift;				// auxiliary ising shift
            ImTG itBins;

            void updateContribution(int sign);
            void updateContribution_OLD(int sign);

            /*! Call Weiss function G0(t1-t2, spin) - alpha(s_ext)_{t1-t2}
             */
            inline RealT g0Call(RealT t1, int spin, int s_ext, RealT t2)
            {
                const RealT t = t1-t2;
                if(t != 0.0) return g0(t,spin);
                return g0(0.0,spin) - (0.5 + (2*(s_ext==spin)-1)*zeroShift);
            }
    };

}
#endif
