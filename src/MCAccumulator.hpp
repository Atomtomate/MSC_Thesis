#ifndef MCACC_HPP_
#define MCACC_HPP_

#include "GreensFct.hpp"
#include "Config.hpp"
#include "IOhelper.hpp"

#include <array>

#include <boost/mpi.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>


#include "fftw3.h"


namespace DMFT
{

    //TODO: implen
    template <int BinSize> //TODO: do a template version of this
        class MCAccumulator
        {
            private:
                using AccT = boost::accumulators::accumulator_set<RealT, boost::accumulators::features<boost::accumulators::tag::sum, boost::accumulators::tag::moment<2> > >;
                int totalSign;
                const GreensFct &g0;
                const Config& config;
                bool dataReady;
                FFT fft;

                std::vector<AccT> binsUp;
                std::vector<AccT> binsDown;

                /*! Saves value at imaginary time tau
                 *
                 *  @param [in]  tau    imaginary time
                 *  @param [in]  val    measurement at time tau
                 *  @param [in]  sign   sign of controbution (>0 for bosonic+WeakCoupling)
                 */
                void push(RealT tau, RealT valDown, RealT valUp, int sign)
                {
                    const int sign2 = 2*(tau>0)-1;
                    int index = static_cast<int>(BinSize * (tau + (tau<0)*config.beta)/config.beta );
                    binsDown[index](sign2*sign*valDown);
                    binsUp[index](sign2*sign*valUp);
                }

                // list of accumulators, one for ech registered process
                // each accumulator has several (exponential?) bins
            public:
                MCAccumulator(const GreensFct &g0, const Config& c):
                    g0(g0), config(c), totalSign(0), dataReady(false), binsUp(BinSize, AccT() ), binsDown(BinSize, AccT() ), fft(c.beta)
            {
            };

                virtual ~MCAccumulator(void) {}

                /*! Collects data from MPI processes.
                 *
                 *  Data is expected to be formatted as std::vector of length n+1.
                 *  The first n/3 elements are imaginary time points,  after that n/3 values for spin DOWN, then spin UP.
                 *  The last element is the sign (+1, -1, 0 to use last sign (sample rejected))
                 */
                void collect(void)
                {
                    while(!config.isGenerator)
                    {
                        boost::mpi::status msg = config.world.probe();
                        if( msg.tag() == static_cast<int>(MPI_MSG_TAGS::DATA))
                        {
                            std::vector<RealT> data;
                            config.world.recv(msg.source(), msg.tag(), data);
                            int dataLength = static_cast<int>((data.size()-1)/3);
                            int sign = data[data.size()-1]; 
                            totalSign += sign;
                            for(int i = 0; i < dataLength; i++)
                                push(data[i],data[dataLength+i],data[2*dataLength+i], sign);
                            /*std::vector<RealT> outVdown(BinSize,0);
                            std::vector<RealT> outVup(BinSize,0);
                            for(int i = 0; i < BinSize; i++)
                            {
                                outVdown[i] = boost::accumulators::sum(binsDown[i]);
                                outVup[i] = boost::accumulators::sum(binsUp[i]);
                                //g0_tmp(i,DOWN) = g0.getByT(t, DOWN);
                                //g0_tmp(i,UP) = g0.getByT(t, UP);
                            }
                            //LOG(INFO) << "-----------------mcacc--------------------: " << totalSign  << "\n"
                            //    << outVdown << "\n" << outVup << "\n"
                            //    << "----------------->mcacc<--------------------";
                            */
                        }
                        else if (msg.tag() == static_cast<int>(MPI_MSG_TAGS::SAMPLING_END))
                        {
                            LOG(INFO) << "sampling end.";
                            //LOG(INFO) << "end comm";
                            config.world.recv(msg.source(), msg.tag());
                            //TODO: sync all collectors
                            dataReady = true;
                            computeImpGF();

                            // cleanup
                            totalSign = 0;
                            dataReady = false;
                            for(int i = 0; i < binsDown.size(); i++)
                            {
                                binsDown[i] = AccT();
                                binsUp[i] = AccT();
                            }
                        }
                        else if(msg.tag() == static_cast<int>(MPI_MSG_TAGS::FINALIZE))
                        {
                            //cleanup?
                            LOG(INFO) << "Closing accumulator";
                            //delete this;
                            ///~MCAccumulator<BinSize>();
                            totalSign = 0;
                            dataReady = false;
                            for(int i = 0; i < binsDown.size(); i++)
                            {
                                binsDown.clear();
                                binsUp.clear();
                            }
                            delete this;
                            return;
                        }
                        else
                        {
                            LOG(ERROR) << "unrecognized message tag: " << msg.tag();
                        }
                    }
                }


                void computeImpGF(void)
                {
                    LOG(INFO) << "Computing new GImp";
                    //TODO: better warning output
                    if(!dataReady)
                        LOG(WARNING) << "data computation before all processes finished";
                    RealT tInc = config.beta/static_cast<RealT>(BinSize);
                    ImTG g0_tmp(BinSize,2);
                    ImTG gImp_tmp_naive(BinSize,2);
                    ImTG gImp_tmp(BinSize,2);
                    ImTG sBins(BinSize,2);

                    for(int i = 0; i < BinSize; i++)
                    {
                        const RealT t = i*tInc;
                        sBins(i,DOWN) = boost::accumulators::sum(binsDown[i]);
                        sBins(i,UP) = boost::accumulators::sum(binsUp[i]);
                        g0_tmp(i,DOWN) = g0(t, DOWN);
                        g0_tmp(i,UP) = g0(t, UP);
                    }
                    std::vector<RealT> mpi_in = fft.conv_naive(g0.getItGF(), sBins, totalSign);
                    boost::mpi::broadcast(config.world, mpi_in, 0);
                    //fft.conv(g0_tmp, sBins, gImp_tmp);
                    //IOhelper::plot(gImp_tmp, config.beta, "gImp convolution");
                    //IOhelper::plot(gImp_tmp_naive, config.beta, "gImp convolution naive");

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
