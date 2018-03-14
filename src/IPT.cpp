#include "IPT.hpp"

namespace DMFT
{
    IPT::IPT(std::string& outDir, GreensFct* const g0, GreensFct* const gImp, const Config& config, const RealT D):
        g0(g0), gImp(gImp), config(config), g0Info( "G0_b"+std::to_string(config.beta)+"_U" + std::to_string(config.U) ), gImpInfo("GImp_b"+std::to_string(config.beta)+"_U" + std::to_string(config.U)), selfE(config.beta, true, true), seLInfo("SelfE_b"+std::to_string(config.beta)+"_U" + std::to_string(config.U)), ioh(outDir, config), D(D)
    {
        ioh.addGF(*g0,   g0Info);
        ioh.addGF(*gImp, gImpInfo);
        ioh.addGF(selfE, seLInfo);

    }

    void IPT::solve(const unsigned int iterations, const unsigned int stat_cycles, const bool out_intermediate, const bool init_insulating)
    {
        if(false)
        {
        for(unsigned int i = 0; init_insulating && (i < 50); i++)
        {
            // set weiss function to semi circular
            for(int f = 0; f < _CONFIG_spins; f++){
                for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
                {
                    int n_g0 = n - (!g0->isSymmetric())*_CONFIG_maxMatsFreq/2;
                    auto iwn = ComplexT(0.,mFreqS(n_g0, config.beta));
                    g0->setByMFreq(n_g0, f, ComplexT(0.,-1./(mFreqS(n_g0, config.beta) - D*D*std::imag(gImp->getByMFreq(n_g0, f))/4.0)));
                }
            }
            g0->markMSet();
            g0->transformMtoT();
            for(int f = 0; f < _CONFIG_spins; f++){
            for(unsigned int it = 0; it < _CONFIG_maxTBins; it++)
            {
                RealT t = config.beta*(it+1)/(_CONFIG_maxTBins+1);
                selfE.setByT(t, f, 20.*(*g0)(t,f)*(*g0)(t,f)*(*g0)(t,f));
            }
            }
            selfE.markTSet();
            selfE.transformTtoM();
            for(int f = 0; f < _CONFIG_spins; f++){
                for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
                {
                    const int n_gi = n - (!gImp->isSymmetric())*_CONFIG_maxMatsFreq/2;
                    gImp->setByMFreq(n_gi, f, ComplexT(0.,1./(1./std::imag(g0->getByMFreq(n_gi ,f)) + std::imag(selfE.getByMFreq(n_gi, f)))));
                }
            }
        }
        }
        for(unsigned int i = 0; i < (iterations + stat_cycles) ; i++)
        {
            for(int f = 0; f < _CONFIG_spins; f++){
                for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
                {
                    int n_g0 = n - (!g0->isSymmetric())*_CONFIG_maxMatsFreq/2;
                    auto iwn = ComplexT(0.,mFreqS(n_g0, config.beta));
                    g0->setByMFreq(n_g0, f, 1./(iwn - D*D*gImp->getByMFreq(n_g0, f)/4.0));
                }
            }
            g0->markMSet();
            //g0->shift(config.U/2.0);
            g0->transformMtoT();

            for(int f = 0; f < _CONFIG_spins; f++){
            for(unsigned int it = 0; it < _CONFIG_maxTBins; it++)
            {
                RealT t = config.beta*(it)/(_CONFIG_maxTBins);
                // shifted Self energy, without config.U/2. for non ph symmetrized
                // config.U/2. -  
                selfE.setByT(t, f, config.U*config.U*(*g0)(t,f)*(*g0)(t,f)*(*g0)(t,f));
            }
            }
            selfE.transformTtoM();

            for(int f = 0; f < _CONFIG_spins; f++){
                for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
                {
                    const int n_gi = n - (!gImp->isSymmetric())*_CONFIG_maxMatsFreq/2;
                    gImp->setByMFreq(n_gi, f, ComplexT( 0. , std::imag(g0->getByMFreq(n_gi ,f)/(ComplexT(1.,0.) - g0->getByMFreq(n_gi ,f)*selfE.getByMFreq(n_gi, f)))));
                }
            }
            if(out_intermediate)
            {
                ioh.setIteration(i+1);
                ioh.writeToFile();
            }
            if(i >= iterations)
            {
                gImp->transformMtoT();
                statAcc(*gImp);
                statAccSE(selfE);
            }
        }
        statAcc.setGF(*gImp);
        statAccSE.setGF(selfE);
        //ioh.writeFinalToFile(*gImp, gImpInfo, false, config.U, false);
        ioh.writeFinalToFile(selfE, seLInfo, true, config.U, false);
    }
}
