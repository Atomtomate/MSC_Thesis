#ifndef STAT_ACC_HPP_
#define STAT_ACC_HPP_

#include "Config.hpp"
#include "GreensFct.hpp"

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
#include <vector>

namespace DMFT
{

template<unsigned int FLAVORS, unsigned int N_MF, unsigned int N_IT>
class StatAcc
{

    using HistAccT = boost::accumulators::accumulator_set< RealT,
          boost::accumulators::features< 
            boost::accumulators::tag::mean,
            boost::accumulators::tag::variance,
            boost::accumulators::tag::skewness,
            boost::accumulators::tag::kurtosis> >;
    public:
        enum class ExtractType : unsigned char {
            MF_REAL = 0,
            MF_IMAG = 1,
            IT_REAL = 2
        };
        enum class StatType : unsigned char {
            MEAN = 0,
            VARIANCE  = 1,
            SKEWNESS = 2,
            KURTOSIS = 3
        };

        void operator()(GreensFct& gf)
        {
            const MatG gwn = gf.g_wn;
            const ImTG git = gf.g_it;
            for(int f = 0; f< FLAVORS; f++)
            {
                for(unsigned int n = 0; n < N_MF; n++)
                {
                    mf_re[n][f](gwn(n,f).real());
                    mf_im[n][f](gwn(n,f).imag());
                }
                for(unsigned int t = 0; t < N_IT; t++)
                {
                    it_re[t][f](git(t,f));
                }
            }
        }

        void setGF(GreensFct& gf)
        {
            MatG gwn = gf.g_wn;
            ImTG git = gf.g_it;
            MatG gwn_std = gf.g_wn;
            ImTG git_std = gf.g_it;
            for(unsigned f = 0; f < FLAVORS; f++)
            {
                auto g_wn_re = extract(ExtractType::MF_REAL, StatType::MEAN, f);
                auto g_wn_im = extract(ExtractType::MF_IMAG, StatType::MEAN, f);
                auto g_wn_var_re = extract(ExtractType::MF_REAL, StatType::VARIANCE, f);
                auto g_wn_var_im = extract(ExtractType::MF_IMAG, StatType::VARIANCE, f);
                auto g_it = extract(ExtractType::IT_REAL, StatType::MEAN, f);
                auto g_it_var = extract(ExtractType::IT_REAL, StatType::VARIANCE, f);
                for(unsigned int n = 0; n < N_MF; n++)
                {
                    gwn(n,f) = ComplexT(g_wn_re[n], g_wn_im[n]);
                    gwn_std(n,f) = ComplexT(std::sqrt(g_wn_var_re[n]), std::sqrt(g_wn_var_im[n]));
                }
                git.col(f) = Eigen::Map<Eigen::VectorXd>(g_it.data(), N_IT);
                git_std.col(f) = Eigen::Map<Eigen::VectorXd>(g_it_var.data(), N_IT);
            }
            git_std = git_std.sqrt();
            gf.g_wn = gwn;
            gf.g_wn_std = gwn_std;
            gf.g_it = git;
            gf.g_it_std = git_std;
        }

        std::vector<RealT> extract(const ExtractType et, const StatType st, const unsigned int f) const
        {
            std::vector<RealT> res;
            unsigned mn = 0;
            switch(et)
            {
                case ExtractType::MF_REAL:
                case ExtractType::MF_IMAG:
                    res.reserve(N_MF);
                    mn = N_MF;
                    break;
                case ExtractType::IT_REAL:
                    res.reserve(N_IT);
                    mn = N_IT;
                    break;
            }
            for(unsigned int n = 0; n < mn; n++)
            {
                RealT val;
                switch(st)
                {
                    case StatType::MEAN:
                        switch(et)
                        {
                            case ExtractType::MF_REAL:
                                val = boost::accumulators::mean(mf_re[n][f]);
                                break;
                            case ExtractType::MF_IMAG:
                                val = boost::accumulators::mean(mf_im[n][f]);
                                break;
                            case ExtractType::IT_REAL:
                                val = boost::accumulators::mean(it_re[n][f]);
                                break;
                        }
                        break;
                    case StatType::VARIANCE:
                        switch(et)
                        {
                            case ExtractType::MF_REAL:
                                val = boost::accumulators::variance(mf_re[n][f]);
                                break;
                            case ExtractType::MF_IMAG:
                                val = boost::accumulators::variance(mf_im[n][f]);
                                break;
                            case ExtractType::IT_REAL:
                                val = boost::accumulators::variance(it_re[n][f]);
                                break;
                        }
                        break;
                    case StatType::SKEWNESS:
                        switch(et)
                        {
                            case ExtractType::MF_REAL:
                                val = boost::accumulators::skewness(mf_re[n][f]);
                                break;
                            case ExtractType::MF_IMAG:
                                val = boost::accumulators::skewness(mf_im[n][f]);
                                break;
                            case ExtractType::IT_REAL:
                                val = boost::accumulators::skewness(it_re[n][f]);
                                break;
                        }
                        break;
                    case StatType::KURTOSIS:
                        switch(et)
                        {
                            case ExtractType::MF_REAL:
                                val = boost::accumulators::kurtosis(mf_re[n][f]);
                                break;
                            case ExtractType::MF_IMAG:
                                val = boost::accumulators::kurtosis(mf_im[n][f]);
                                break;
                            case ExtractType::IT_REAL:
                                val = boost::accumulators::kurtosis(it_re[n][f]);
                                break;
                        }
                        break;
                }
                res.push_back(val);
            }
            return res;
        }


    private:
        std::array<std::array<HistAccT,FLAVORS>, N_MF> mf_im;
        std::array<std::array<HistAccT,FLAVORS>, N_MF> mf_re;
        std::array<std::array<HistAccT,FLAVORS>, N_IT> it_re;
};
}
#endif
