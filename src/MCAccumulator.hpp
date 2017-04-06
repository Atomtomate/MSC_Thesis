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
                int totalSign;
                int lastSign;
                const GreensFct &g0;
                const Config& config;
                bool dataReady;

                AccT binsUp;
                std::vector<AccT> binsDown;

                /*! Saves value at imaginary time tau
                 *
                 *  @param [in]  tau    imaginary time
                 *  @param [in]  val    measurement at time tau
                 *  @param [in]  sign   sign of controbution (>0 for bosonic+WeakCoupling)
                 */
                void push(RealT tau, RealT valUp, RealT valDown, int sign)
                {
                    if(!sign) sign = lastSign;		    // update got rejected, use last sign
                    lastSign = sign;			    // remember last sign
                    totalSign += 1;//sign;
                    const int sign2 = 2*(tau>0)-1;
                    int index = static_cast<int>(BinSize * (tau + (tau<0)*config.beta) );
                    binsUp(sign2*sign*valUp);
                    //binsDown[0](sign2*sign*valDown);
                }

                // list of accumulators, one for ech registered process
                // each accumulator has several (exponential?) bins
            public:
                MCAccumulator(const GreensFct &g0, const Config& c):
                    g0(g0), config(c), totalSign(0), lastSign(0), dataReady(false)
                                    {
                };

                /*! Collects data from MPI processes.
                 *
                 *  Data is expected to be formatted as std::vector of length n+1.
                 *  The first n/3 elements are imaginary time points,  after that n/3 values for spin UP, then spin DOWN.
                 *  The last element is the sign (+1, -1, 0 to use last sign (sample rejected))
                 */
                void collect(void)
                {
                    while(!config.isGenerator)
                    {
                        boost::mpi::status msg = config.world.probe();
                        if( msg.tag() == static_cast<int>(MPI_MSG_TAGS::DATA))
                        {
                            boost::optional<std::vector<RealT> > data;
                            data = std::vector<RealT>();
                            config.world.recv(msg.source(), msg.tag(), *data);
                            do
                            {
                                if(data)
                                {
                                    std::vector<RealT> vec = *data;
                                    int dataLength = static_cast<int>((vec.size()-1)/3);
                                    int sign = vec[3*dataLength]; 
                                    for(int i = 0; i < dataLength; i++)
                                        push(vec[i],vec[dataLength+i],vec[2*dataLength+i], sign);
                                    if(dataLength == 0)
                                    {
                                        if(!sign) sign = lastSign;
                                        lastSign = sign;
                                        totalSign += 1;
                                    }
                                }
                            } while(data);
                        }
                        else if (msg.tag() == static_cast<int>(MPI_MSG_TAGS::COMM_END))
                        {
                            config.world.recv(msg.source(), msg.tag());
                            //TODO: sync all collectors
                            dataReady = true;
                            break;
                        }
                    }
                }


                void computeImpGF(void)
                {
                    //TODO: better warning output
                    if(!dataReady)
                        LOG(WARNING) << "data computation before all processes finished";
                    RealT tInc = config.beta/static_cast<RealT>(config.itCount);
                    for(int i = 0; i < config.itCount; i++)
                    {
                        RealT tau = tInc*i;
                        
                    }
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
