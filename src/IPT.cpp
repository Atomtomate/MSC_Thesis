#include "IPT.hpp"

namespace DMFT
{
    IPT::IPT(std::string& outDir, GreensFct* const g0, GreensFct* const gImp, const Config& config, const RealT D):
        g0(g0), gImp(gImp), config(config), g0Info( "G0" + std::to_string(config.U) ), gImpInfo("GImp_U" + std::to_string(config.U)), selfE(config.beta, true, false), seLInfo("SelfE_U" + std::to_string(config.U)), ioh(outDir, config), D(D)
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
                    int n_g0 = n - (!g0->isSymmetric())*_CONFIG_maxMatsFreq/2;
                    auto iwn = ComplexT(0.,mFreqS(n_g0, config.beta));
                    g0->setByMFreq(n_g0, f, ComplexT(0.,-1./(mFreqS(n_g0, config.beta) - D*D*std::imag(gImp->getByMFreq(n_g0, f))/4.0)));
                }
            }
            g0->markMSet();
            //g0->shift(config.U/2.0);
            g0->transformMtoT();

            for(int f = 0; f < _CONFIG_spins; f++){
            for(unsigned int it = 0; it < _CONFIG_maxTBins; it++)
            {
                RealT t = config.beta*(it+1)/(_CONFIG_maxTBins+1);
                // shifted Self energy, without config.U/2. for non ph symmetrized
                // config.U/2. -  
                selfE.setByT(t, f, config.U*config.U*(*g0)(t,f)*(*g0)(t,f)*(*g0)(t,f));
            }
            }
            selfE.markTSet();
            selfE.transformTtoM();
            if(out_intermediate)
            {
                ioh.setIteration(i+1);
                ioh.writeToFile();
            }
            exit(0);

            for(int f = 0; f < _CONFIG_spins; f++){
            for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
            {
                const int n_gi = n - (!gImp->isSymmetric())*_CONFIG_maxMatsFreq/2;
                gImp->setByMFreq(n_gi, f, ComplexT(0.,1./(1./std::imag(g0->getByMFreq(n_gi ,f)) + std::imag(selfE.getByMFreq(n_gi, f)))));
            }
            }
        }
        ioh.writeFinalToFile(*gImp, gImpInfo);
        ioh.writeFinalToFile(selfE, seLInfo, true, config.U);
            gImp->transformMtoT();
            IOhelper::plot(selfE, config.beta, "Self energy");
            IOhelper::plot(*g0, config.beta, "Weiss Function");
            IOhelper::plot(*gImp, config.beta, "Imp Function");
            exit(0);
    }
}
