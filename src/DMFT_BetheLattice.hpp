#ifndef DMFT_BETHE_HPP_
#define DMFT_BETHE_HPP_
#include "GreensFct.hpp"
#include "WeakCoupling.hpp"
#include "StrongCoupling.hpp"
#include "KInt.hpp"
#include "IOhelper.hpp"
#include "ImpSolver.hpp"
#include "integrals.hpp"
#include "StatAcc.hpp"

namespace DMFT
{

template<class ImpSolver>
class DMFT_BetheLattice
{
    //"GImp_b"+std::to_string(config.beta)+"_U" + std::to_string(config.U)   _b"+std::to_string(config.beta)+"_U" + std::to_string(config.U)
    public:
        DMFT_BetheLattice(std::string& outDir, const Config config, RealT mixing, ImpSolver& solver, GreensFct* G0, GreensFct * GImp, const RealT D, bool useHyb = false):
            config(config), mixing(mixing), iSolver(solver), useHyb(useHyb), \
            g0(G0), g0Info( useHyb ? "Hyb" : "G0" ), gImp(GImp), gImpInfo("GImp"), selfE(config.beta, true, true), seLInfo("SelfE"),\
            ioh(outDir, config), D(D)
        {
            std::string tmp( useHyb ? "Hyb_Guess" : "G0_Guess");
            ioh.writeToFile(*g0,tmp);
            ioh.addGF(*g0,   g0Info);
            ioh.addGF(*gImp, gImpInfo);
            ioh.addGF(selfE, seLInfo);
        }

        void solve(const unsigned int iterations, const unsigned long long updates, bool symmetricG0 = false, bool writeOut = true)
        {
            if(config.isGenerator)
            {
                //IOhelper::plot(*g0, config.beta, "Weiss Function");
                //TODO: move IO from generator to accumulators
                for(unsigned int dmftIt = 1;dmftIt < iterations+1; dmftIt++)
                {  
                    for(long unsigned int i=1; i <= 5; i++){
                        iSolver.update(updates/5.0);
                        LOG(INFO) << "MC Walker [" << config.world.rank() << "] at "<< " (" << (20*i) << "%) of iteration " << dmftIt << ". expansion order: " << iSolver.expansionOrder();
                    }
                    ioh.setIteration(dmftIt);
                    if(config.local.rank() == 0)
                    {
                        LOG(INFO) << "finished sampling. average expansion order: " << iSolver.expansionOrder();
                        LOG(INFO) << "measuring impurity Greens function";
                    }
                    iSolver.computeImpGF();
                    if(config.local.rank() == 0)
                    {
                        LOG(INFO) << "forcing paramagnetic solution";
                        LOG(INFO) << "Writing results";
                    }
                    if(!useHyb)
                    {
                        g0->shift(config.U/2.0);
                        g0->markMSet();
                        g0->transformMtoT();
                    }
                    if(config.local.rank() == 0) LOG(INFO) << "Computing new Weiss Green's function";
                    //TODO: vectorize
                    MatG tmpG0(_CONFIG_maxTBins,2);
                    MatG sImp(_CONFIG_maxMatsFreq, _CONFIG_spins);
                    ImTG sImp_it(_CONFIG_maxTBins, _CONFIG_spins);
                    IOhelper::plot(*gImp, config.beta, "Imp GF iteration " + std::to_string(dmftIt));
                    exit(0);
                    // DMFT equation
                    for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                        const int n_g0 = n + ((int)g0->isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                        const int n_se = n + ((int)selfE.isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                        const ComplexT iwn_se(0., mFreqS(n_se, config.beta));
                        for(int s=0;s<_CONFIG_spins;s++){
                            if(useHyb)
                            {
                                /*if(s == 0){
                                LOG(ERROR) << n << ": " << iwn_se << " - " <<  (D*D)*gImp->getByMFreq(n_se, s)/(4.0)<< " - " << 1.0/gImp->getByMFreq(n_se, s) << \
                                    " = " << iwn_se -  (D*D)*gImp->getByMFreq(n_se, s)/(4.0) - 1.0/gImp->getByMFreq(n_se, s);
                                LOG(ERROR) << "AL: " << mFreqS(n_se, config.beta) << " - " << ((D*D)*gImp->getByMFreq(n_se, s)/(4.0)).imag() << " + " << 1.0/((gImp->getByMFreq(n_se, s)).imag()) \
                                    << " = " << ComplexT(0., mFreqS(n_se, config.beta) - ((D*D)*gImp->getByMFreq(n_se, s)/(4.0)).imag() + 1.0/((gImp->getByMFreq(n_se, s)).imag()));
                                }
                                if(n > 20)
                                    exit(0);
                                    */
                                selfE.setByMFreq(n_se, s, iwn_se -  (D*D)*gImp->getByMFreq(n_se, s)/(4.0) - 1.0/gImp->getByMFreq(n_se, s) );
                            }
                            ComplexT tmp = ComplexT(0., mFreqS(n_g0, config.beta)) - (D/2.0)*(D/2.0)*gImp->getByMFreq(n_g0, s);
                            VLOG(5) << n << "=> "<< n_g0 << ": " << ComplexT(config.mu, mFreqS(n_g0,config.beta)) << " - " << (D/2.0)*(D/2.0)*gImp->getByMFreq(n_g0,s) << " = " << tmp;
                            g0->setByMFreq(n_g0, s, 1.0/tmp );
                        }
                    }
                    g0->setParaMagnetic();
                    if(useHyb)
                    {
                        LOG(INFO) << "tf g0";
                        g0->markMSet();
                        g0->transformMtoT();
                        selfE.markMSet();
                        selfE.transformMtoT();
                    }
                    else
                    {
                        //only for WeakCoupling., this is now in update itself
                        sImp = (g0->getMGF().cwiseInverse() - gImp->getMGF().cwiseInverse());
                        selfE.setByMFreq(sImp);
                        selfE.markMSet();
                        selfE.transformMtoT();
                        //g0->shift(config.U/2.0);
                        //g0->markMSet();
                        //g0->transformMtoT();
                    }
                    if(dmftIt == 1)
                    {
                        ioh.setIteration(0);
                        //ioh.writeToFile();
                    }
                    
                    IOhelper::plot(*g0, config.beta, "Hyb Fct iteration " + std::to_string(dmftIt));
                    IOhelper::plot(selfE, config.beta, "Self Energy iteration " + std::to_string(dmftIt));
                    exit(0);
                    //
                    if(writeOut)
                    {
                        ioh.writeToFile();
                    }
                    if(dmftIt == iterations)
                    {
                        statAcc(*gImp);
                        statAccSE(selfE);
                    }

                }
                iSolver.avgN().writeResults(gImpInfo.filename);
                statAcc.setGF(*gImp);
                statAccSE.setGF(selfE);
                ioh.writeFinalToFile(*gImp, gImpInfo, false, config.U, false);
                ioh.writeFinalToFile(selfE, seLInfo, true, config.U, false);

                //config.world.send(0, static_cast<int>(MPI_MSG_TAGS::FINALIZE), 1 );
            }
            else
            {
                //MCAccumulator<_CONFIG_maxSBins> *mcAcc = new MCAccumulator<_CONFIG_maxSBins>(*g0,config);
                //mcAcc->collect();
                //TODO: compute GImp
                //TODO: compute new g0
                //TODO: distribute g0 to generators
            }


        }

        void setDir(std::string d) { ioh.initDir(d); }

    private:
        // general settings
        const Config	config;
        IOhelper	    ioh;

        // lattice specific
        // TODO separate class
        RealT D;
        const bool useHyb;
        RealT mixing;

        ImpSolver& iSolver;
        LogInfos	g0Info;
        LogInfos	gImpInfo;
        LogInfos	seLInfo;
        GreensFct*  g0;
        GreensFct*	gImp;
        GreensFct   selfE;
        StatAcc<_CONFIG_spins, _CONFIG_maxMatsFreq, _CONFIG_maxTBins> statAcc;
        StatAcc<_CONFIG_spins, _CONFIG_maxMatsFreq, _CONFIG_maxTBins> statAccSE;
};
} // end namespace DMFT
#endif
