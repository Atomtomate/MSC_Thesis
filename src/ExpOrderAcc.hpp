#ifndef EXP_ORD_ACC_HPP_
#define EXP_ORD_ACC_HPP_

#include "Config.hpp"
#include "IOhelper.hpp"

#include <boost/serialization/vector.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/skewness.hpp>

#include <iostream>
#include <tuple>
#include <algorithm>
#include <sstream> 
#include <string>
namespace DMFT
{

template<unsigned int FLAVORS>
class ExpOrderAcc
{
    using HistAccT = boost::accumulators::accumulator_set< RealT,
          boost::accumulators::features< 
            boost::accumulators::tag::mean,
            boost::accumulators::tag::variance,
            boost::accumulators::tag::skewness,
            boost::accumulators::tag::kurtosis,
            boost::accumulators::tag::density> >;
    typedef boost::iterator_range<std::vector<std::pair<RealT, RealT> >::iterator > HistogramT;

    public:
        ExpOrderAcc(const Config *const c): c(c), ioh(*c), expansionOrderAcc(FLAVORS, HistAccT(boost::accumulators::tag::density::num_bins = 20, boost::accumulators::tag::density::cache_size = 10000)) 
        {}
        
        inline void operator()(RealT val, unsigned int f = 0)
        {
            expansionOrderAcc[f](val);
        }

        void writeResults(std::string name = "")
        {
            std::array<HistogramT, FLAVORS> histL;
            RealT mean = 0, var = 0, skew = 0, kurt = 0;
            for(int f = 0; f < FLAVORS; f++)
            {
                histL[f] = boost::accumulators::density(expansionOrderAcc[f]);
                mean += boost::accumulators::mean(expansionOrderAcc[f]);
                var += boost::accumulators::variance(expansionOrderAcc[f]);
                skew += boost::accumulators::skewness(expansionOrderAcc[f]);
                kurt += boost::accumulators::kurtosis(expansionOrderAcc[f]);
            }
            //TODO: this is a workaround, assuming <n_up> ~ <n_down> and
            //FLAVORS == 2
            if(FLAVORS > 2) LOG(WARNING) << "histogram computation for more than two flavors not yet implemented";
            var /= FLAVORS; skew /= FLAVORS; kurt /= FLAVORS;

            int i1 = 0, i2 = 0;
            int i = 0;
            std::vector<std::pair<RealT, RealT>> histOut;
            while(i < histL[0].size())
            {
                    histOut.push_back(std::make_pair(histL[0][i1].first, histL[0][i1].second));
                    i1 += 1; i+=1;
            }
            std::stringstream ss;
            std::string filename = name + "ExpansionOrder_U"+std::to_string(c->U);
            ss << "# mean\tvariance\tskewness\tkurtosis" << std::endl << mean << "\t" << var << "\t" << skew << "\t" << kurt << std::endl; 
            ss << "# l bound\tcount" << std::endl;
            for(int i =0; i < histOut.size(); i++)
            {
                ss << histOut[i].first << "\t" << histOut[i].second << std::endl;
            }
            ioh.writeToFile(ss.str(), filename);
        }

    private:
        const Config * const c;
        IOhelper ioh;
        std::vector<HistAccT> expansionOrderAcc;

};


}   // end namspace DMFT
#endif
