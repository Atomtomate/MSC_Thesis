#ifndef MCACC_HPP_
#define MCACC_HPP_



#include <array>

#include <boost/mpi.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

#include "Config.hpp"

namespace DMFT
{

template <int NBins> //TODO: do a template version of this
class MCAccumulator
{
    private:
        using AccT = boost::accumulators::accumulator_set<RealT, boost::accumulators::tag::mean, boost::accumulators::tag::moment<2> >;
        std::array<AccT, NBins> bins;
        Config& config;
        // list of accumulators, one for ech registered process
        // each accumulator has several (exponential?) bins
    public:
        MCAccumulator(Config& c): config(c) {};

        /*! Saves value at imaginary time tau
         *
         *  @param [in]  tau    imaginary time
         *  @param [in]  val    measurement at time tau
         */
        void push(RealT tau, RealT val)
        {
            bins[static_cast<int>(_CONFIG_maxSBins * (tau + (tau<0)*config.beta)/config.beta )](val);
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
