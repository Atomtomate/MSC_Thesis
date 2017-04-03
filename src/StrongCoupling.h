#ifndef STRONG_COUPLING_H_
#define STRONG_COUPLING_H_


#include "Config.hpp"
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
	typedef std::pair<const double,const char> SConfig;
	struct compare {
  		inline bool operator()(const SConfig& lhs, const SConfig& rhs) const {
          	return (lhs.first) < (rhs.first);}
	};
	typedef std::set<SConfig, compare> SConfigTOL;                                  // Time orderer list (set)
	typedef std::vector<SConfig> SConfigSOL;                                        // Ordered by time of generation

  #ifdef TIME_ORDERED_CONFIG_LIST
        typedef SConfigTOL SConfigL;                                            // TODO: switch to Boost::set
        inline static void pushConfig(SConfigL& l, const SConfig& c){
                l.insert(c);
        }
        inline static void deleteConfig(SConfigL& l, const int pos){
                auto sConfIt = confs.begin();                                   // TODO: make this a class var?
                std::advance(sConfIt, pos);
                auto tmpConf = *sConfIt;
                confs.erase(sConfIt);
        }

        typedef
  #else                                                                           // TODO: reseve size
        typedef SConfigSOL SConfigL;
        inline static void pushConfig(SConfigL& l, const SConfig& c){
                l.push_back(c);
        }
        inline static void deleteConfig(SConfigL& l, const int pos){
                l.pop_back();
        }
  #endif
        public:
                StrongCoupling(SConfigL& confs, MatrixT& hybr, const RealT U, const RealT zeroShift, const RealT beta, const unsigned int burninSteps);
                ~StrongCoupling();


                /**     Updates the time ordered spin configuration and the inverse
                 *  weiss greens function M by inserting or removing one configuration.
                 *  @param C List of configurations
                 *  @param M Inverse Weiss Green's function
                 *  @param G0 Weiss Green's function
                 *  @param U Interaction
                 *  @param beta Inverse temperature
                 */
                int update(const RealT U,const RealT beta);
                void doMeasurement(void);

                /** computes the acceptance rate by evaluation of the determinant
                 * ratios \f$ \frac{-\beta U}{(n+1}  \prod \limits_{\sigma}
                 * \frac{det[M(n+1,\sigma)^{-1}]}{det[M(n,\sigma)^{-1}]}\f \quad GKW (8.36)$
                 */
                RealT acceptanceR(const RealT U, const RealT beta) const;

        private:
                MatrixT& hybr;
                SConfigL& confs;
                unsigned int steps;                     // number of updates
                const unsigned int burninSteps;
                unsigned int n;                         // expansion order (number of used rows/cols)
                const RealT zeroShift;                  // auxiliary ising shift
                //Hybr g0;                                  // Hybridization function
};

} //end namespace DMFT
#endif

