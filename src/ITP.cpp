#include "IPT.hpp"

namespace DMFT
{
    IPT::IPT(GreensFct &g0, GreensFct &gImp, const Config& config): g0(g0), gImp(gImp), config(config) {}
    void IPT::update(void)
    {
        auto tmp = config.U * config.U * g0.getItGF() * g0.getItGF() * g0.getItGF(); 
    }
}
