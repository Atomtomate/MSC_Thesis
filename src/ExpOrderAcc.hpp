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
            boost::accumulators::tag::count,
            boost::accumulators::tag::mean,
            boost::accumulators::tag::variance,
            boost::accumulators::tag::skewness,
            boost::accumulators::tag::kurtosis,
            boost::accumulators::tag::density> >;
    typedef boost::iterator_range<std::vector<std::pair<RealT, RealT> >::iterator > HistogramT;

    public:
        ExpOrderAcc(const Config c): c(c), expansionOrderAcc(FLAVORS, HistAccT(boost::accumulators::tag::density::num_bins = 20, boost::accumulators::tag::density::cache_size = 10000)) 
        {}

        virtual ~ExpOrderAcc(){
            //LOG(ERROR) << "destroying exp order acc";
        }
        
        inline void operator()(RealT val, unsigned int f = 0)
        {
            //LOG(ERROR) << boost::accumulators::count(expansionOrderAcc[f]) << ", " << val << " : " << f;
            expansionOrderAcc[f](val);
        }

        std::vector<RealT> get_stat(const int f)
        {
            if(f >= FLAVORS) LOG(ERROR) << "Invalid flavor!";
            std::vector<RealT> res(4);
            LOG(ERROR) << "deactivated";
            return res;
            res[0] = boost::accumulators::mean(expansionOrderAcc[f]);
            res[1] = boost::accumulators::variance(expansionOrderAcc[f]);
            res[2] = boost::accumulators::skewness(expansionOrderAcc[f]);
            res[3] = boost::accumulators::kurtosis(expansionOrderAcc[f]);
            return res;
        }

        void writeResults(IOhelper ioh, std::string name = "")
        {
            std::stringstream ss;
            std::string filename = name + "ExpansionOrderb_"+std::to_string(c.beta)+"_U"+std::to_string(c.U);
            for(int f = 0; f < FLAVORS; f++)
            {
                HistogramT histL = boost::accumulators::density(expansionOrderAcc[f]);
                RealT mean = boost::accumulators::mean(expansionOrderAcc[f]);
                RealT var = boost::accumulators::variance(expansionOrderAcc[f]);
                RealT skew = boost::accumulators::skewness(expansionOrderAcc[f]);
                RealT kurt = boost::accumulators::kurtosis(expansionOrderAcc[f]);
                std::vector<std::pair<RealT, RealT>> histOut;
                ss << "# Flavor: " << f << std::endl;
                ss << "# mean\tvariance\tskewness\tkurtosis" << std::endl << mean << "\t" << var << "\t" << skew << "\t" << kurt << std::endl; 
                ss << "# l bound\tcount" << std::endl;
                for(int i =0; i < histL.size(); i++)
                {
                    ss << histL[i].first << "\t" << histL[i].second << std::endl;
                }
            }
            ioh.writeToFile(ss.str(), filename);
        }

        void reset()
        {
           // LOG(WARNING) << "resetting expansion order accumulator";
            for(int f = 0; f < FLAVORS; f++)
            {
                expansionOrderAcc[f] = HistAccT(boost::accumulators::tag::density::num_bins = 20, boost::accumulators::tag::density::cache_size = 10000);
            }
        }

    private:
        const Config c;
        std::vector<HistAccT> expansionOrderAcc;

};


}   // end namspace DMFT
#endif
