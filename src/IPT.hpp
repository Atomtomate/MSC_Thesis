#ifndef IPT_HPP_
#define IPT_HPP_

#include "Config.hpp"
#include "GreensFct.hpp"
#include "ImpSolver.hpp"
#include "IOhelper.hpp"

namespace DMFT
{
class IPT
{
    public:
        IPT(std::string& outDir, GreensFct* const g0, GreensFct* const gImp, const Config& config, const RealT D);
        void solve(unsigned int iterationsi = 10, bool out_intermediate = false);

    private:
        const Config& config;
        IOhelper	    ioh;

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
