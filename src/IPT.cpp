#include "IPT.hpp"

namespace DMFT
{
    IPT::IPT(std::string& outDir, GreensFct* const g0, GreensFct* const gImp, const Config& config, const RealT D):
        g0(g0), gImp(gImp), config(config), g0Info( "G0" ), gImpInfo("GImp"), selfE(config.beta, true, false), seLInfo("SelfE"), ioh(outDir, config), D(D)
    {
        ioh.addGF(*g0,   g0Info);
        ioh.addGF(*gImp, gImpInfo);
        ioh.addGF(selfE, seLInfo);

    }

    void IPT::solve(unsigned int iterations, bool out_intermediate)
    {
        for(unsigned int i = 0; i < iterations; i++)
        {
            for(int f = 0; f < _CONFIG_spins; f++){
            for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
            {
                const int n_g0 = n + ((int)g0->isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                const ComplexT iwn_g0(0., mFreqS(n_g0, config.beta));
                g0->setByMFreq(n_g0, f, 1./(iwn_g0 - D*D*gImp->getByMFreq(n_g0, f)/4.0));
            }
            }
            g0->markMSet();
            g0->shift(config.U/2.0);
            g0->transformMtoT();

            for(int f = 0; f < _CONFIG_spins; f++){
            for(unsigned int it = 0; it < _CONFIG_maxTBins; it++)
            {
                RealT t = config.beta*it/_CONFIG_maxTBins;
                // shifted Self energy, without config.U/2. for non ph symmetrized
                // config.U/2. -  
                selfE.setByT(t, f, config.U*config.U*(*g0)(t,f)*(*g0)(t,f)*(*g0)(t,f));
            }
            }
            selfE.markTSet();
            selfE.transformTtoM();

            for(int f = 0; f < _CONFIG_spins; f++){
            for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
            {
                const int n_gi = n + ((int)gImp->isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                const ComplexT iwn_gi(0., mFreqS(n_gi, config.beta));
                gImp->setByMFreq(n_gi, f, 1./(1./g0->getByMFreq(n_gi ,f) - selfE.getByMFreq(n_gi, f)));
            }
            }
            if(out_intermediate)
            {
                ioh.setIteration(i+1);
                ioh.writeToFile();
            }
        }
        ioh.writeFinalToFile(*gImp, gImpInfo);
    }
}
