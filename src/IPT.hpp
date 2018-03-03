#ifndef IPT_HPP_
#define IPT_HPP_

#include "Config.hpp"
#include "GreensFct.hpp"
#include "ImpSolver.hpp"
#include "IOhelper.hpp"
#include "StatAcc.hpp"

namespace DMFT
{
class IPT
{
    public:
        IPT(std::string& outDir, GreensFct* const g0, GreensFct* const gImp, const Config& config, const RealT D);
        void solve(const unsigned int iterations = 10, const unsigned int stat_cycles = 100, const bool out_intermediate = false, const bool init_insulating = false);

    private:
        const Config& config;
        IOhelper	    ioh;
        StatAcc<_CONFIG_spins, _CONFIG_maxMatsFreq, _CONFIG_maxTBins> statAcc;
        StatAcc<_CONFIG_spins, _CONFIG_maxMatsFreq, _CONFIG_maxTBins> statAccSE;

        LogInfos	g0Info;
        LogInfos	gImpInfo;
        LogInfos	seLInfo;

        GreensFct* const g0;
        GreensFct* const gImp;
        GreensFct selfE;

        const RealT D;
};

}   // end namespace DMFT
#endif
