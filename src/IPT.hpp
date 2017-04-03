#ifndef IPT_HPP_
#define IPT_HPP_

#include "GreensFct.hpp"
#include "Config.hpp"

namespace DMFT
{
class IPT
{
    public:
        IPT(GreensFct &g0, GreensFct &gImp, const Config& config);
        void update(void);

    private:
        GreensFct& g0;
        GreensFct& gImp;
        const Config& config;
};

}   // end namespace DMFT
#endif
