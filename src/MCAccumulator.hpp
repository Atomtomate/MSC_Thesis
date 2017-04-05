#ifndef MCACC_HPP_
#define MCACC_HPP_

#include "GreensFct.hpp"

#include <array>

#include <boost/mpi.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

#include "Config.hpp"

namespace DMFT
{

    template <int BinSize> //TODO: do a template version of this
        class MCAccumulator
        {
            private:
                using AccT = boost::accumulators::accumulator_set<RealT, boost::accumulators::tag::mean, boost::accumulators::tag::moment<2> >;
                AccT bins;
                int totalSign;
                int lastSign;
                const GreensFct &g0;
                const Config& config;

                // list of accumulators, one for ech registered process
                // each accumulator has several (exponential?) bins
            public:
                MCAccumulator(const GreensFct &g0, const Config& c): g0(g0), config(c), totalSign(0), lastSign(0)
                {
                };

                void collect(void)
                {
                    while(!config.isGenerator)
                    {
                        boost::mpi::status msg = config.world.probe();
                        if( msg.tag() == static_cast<int>(MPI_MSG_TAGS::DATA))
                        {
                            
                        }
                        else if (msg.tag() == static_cast<int>(MPI_MSG_TAGS::COMM_END))
                        {
                            config.world.recv(msg.source(), msg.tag());
                            //TODO: sync all collectors
                        }
                    }
                }

                /*! Saves value at imaginary time tau
                 *
                 *  @param [in]  tau    imaginary time
                 *  @param [in]  val    measurement at time tau
                 *  #param [in]  sign   sign of controbution (>0 for bosonic+WeakCoupling)
                 */
                void push(RealT tau, RealT val, int sign)
                {
                    if(!sign) sign = lastSign;		    // update got rejected, use last sign
                    lastSign = sign;			    // remember last sign
                    totalSign += 1;//sign;
                    const int sign2 = 2*(tau>0)-1;
                    //[static_cast<int>(BinSize * (tau + (tau<0)*config.beta) )]
                    bins(sign2*sign*val);
                }

                void computeImpGF(boost::mpi::communicator &c)
                {
                    for(int n=0;n<_CONFIG_maxMatsFreq; n++){
                        /*   RealT t 		= config.beta*static_cast<RealT>(n)/_CONFIG_maxMatsFreq;
                             ComplexT sumWn	= 0.0;
                             RealT sumIt		= 0.0;
                             RealT mfreq		= mFreq(n,config.beta);
                             for(int j=0; j<_CONFIG_maxSBins; j++){
                             const RealT bVal = bins(j);
                             if(bVal == 0.0) continue;
                             const RealT tp = static_cast<RealT>(j)/_CONFIG_maxSBins;
                             sumWn += std::exp(ComplexT(0.0, mfreq*tp))*bVal;
                             sumIt += g0(t-tp,s)*bVal;
                             }
                        //TODO: better accumulator
                        //if(std::isnan(std::real(sumWn)) or std::isinf(std::real(sumWn)) or std::isnan(std::imag(sumWn)) or std::isinf(std::imag(sumWn)) ) LOG(ERROR) << "Overflow during computation of G_Imp(i wn)";
                        //if(std::isnan(sumIt) or std::isinf(sumIt)) LOG(ERROR) << "Overflow during computation of G_Imp(tau)";
                        gImp.setByMFreq(n,s, g0.getByMFreq(n,s) - g0.getByMFreq(n,s)*sumWn/static_cast<RealT>(totalSign));
                        gImp.setByT(t,s, g0(t,s) - sumIt/totalSign); */
                    }
                }

                // register thread - generate accumulator
                // get output (generate statistics)

        };

} //end namespace DMFT
#endif
