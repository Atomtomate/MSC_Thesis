#ifndef DMFT_BETHE_HPP_
#define DMFT_BETHE_HPP_
#include "GreensFct.hpp"
#include "WeakCoupling.hpp"
#include "StrongCoupling.hpp"
#include "KInt.hpp"
#include "IOhelper.hpp"
#include "ImpSolver.hpp"
#include "integrals.hpp"

namespace DMFT
{

    template<class ImpSolver>
        class DMFT_BetheLattice
        {
            public:
                DMFT_BetheLattice(std::string& outDir, const Config config, RealT mixing, ImpSolver& solver, GreensFct* G0, GreensFct * GImp, const RealT D, bool useHyb = false):
                    config(config), mixing(mixing), iSolver(solver), useHyb(useHyb), \
                    g0(G0), g0Info( useHyb ? "Hyb" : "G0" ), gImp(GImp), gImpInfo("GImp"), selfE(config.beta, true, false), seLInfo("SelfE"),\
                    ioh(outDir, config), D(D)
                {
                    std::string tmp( useHyb ? "Hyb_Guess" : "G0_Guess");
                    ioh.writeToFile(*g0,tmp);
                    ioh.addGF(*g0,   g0Info);
                    ioh.addGF(*gImp, gImpInfo);
                    ioh.addGF(selfE, seLInfo);
                }

                void solve(const unsigned int iterations, const unsigned long long updates, bool symmetricG0 = false)
                {
                    if(config.isGenerator)
                    {
                        //IOhelper::plot(*g0, config.beta, "Weiss Function");
                        //TODO: move IO from generator to accumulators
                        for(unsigned int dmftIt = 1;dmftIt < iterations+1; dmftIt++)
                        {  
                            if(config.local.rank() == 0) LOG(INFO) << "Computing new Weiss Green's function";
                            //TODO: vectorize
                            MatG tmpG0(_CONFIG_maxTBins,2);
                            MatG sImp(_CONFIG_maxMatsFreq, _CONFIG_spins);
                            ImTG sImp_it(_CONFIG_maxTBins, _CONFIG_spins);
                            // DMFT equation
                            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                                const int n_g0 = n + ((int)g0->isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                                const int n_se = n + ((int)selfE.isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                                const ComplexT iwn_se(0., mFreqS(n_se, config.beta));
                                for(int s=0;s<_CONFIG_spins;s++){
                                    // set SC condition, enforcing Paramagnetic solution
                                    // use symmetry here
                                    if(useHyb)
                                    {
                                        g0->setByMFreq(n_g0, s, + 4.*gImp->getByMFreq(n_g0,s)/(D*D));
                                        //g0->setByMFreq(n_g0, s, - 4.0*gImp->getByMFreq(n_g0,s)/(D*D));
                                        selfE.setByMFreq(n_se, s, iwn_se - 1.0/gImp->getByMFreq(n_se, s) );
                                    }
                                    else
                                    {
                                        ComplexT tmp = ComplexT(config.mu, mFreqS(n_g0,config.beta)) - (D/2.0)*(D/2.0)*gImp->getByMFreq(n_g0,s);
                                        VLOG(5) << n << "=> "<< n_g0 << ": " << ComplexT(config.mu, mFreqS(n_g0,config.beta)) << " - " << (D/2.0)*(D/2.0)*gImp->getByMFreq(n_g0,s) << " = " << tmp;
                                        g0->setByMFreq(n_g0,s, 1.0/tmp );
                                    }
                                }
                            }

                            //IOhelper::plot(*g0, config.beta, "Hyb Fct iteration " + std::to_string(dmftIt));
                            //IOhelper::plot(*gImp, config.beta, "Imp GF iteration " + std::to_string(dmftIt));
                            //exit(0);

                            g0->setParaMagnetic();
                            if(!useHyb)
                            {
                                //only for WeakCoupling., this is now in update itself
                                sImp = (g0->getMGF().cwiseInverse() - gImp->getMGF().cwiseInverse());
                                selfE.setByMFreq(sImp);
                                selfE.setParaMagnetic();
                                g0->shift(config.U/2.0);
                                g0->markMSet();
                                g0->transformMtoT();
                            }
                            selfE.markMSet();
                            selfE.transformMtoT();
                            if(dmftIt == 1)
                            {
                                ioh.setIteration(0);
                                ioh.writeToFile();
                            }
                            
                            //IOhelper::plot(gImp, config.beta, "before measure gImp Function " + std::to_string(dmftIt));

                            //TODO improve this naive loop
                            for(long unsigned int i=1; i <= 20; i++){
                                iSolver.update(updates/20.0);
                                LOG(INFO) << "MC Walker [" << config.world.rank() << "] at "<< " (" << (5*i) << "%) of iteration " << dmftIt << ". expansion order: " << iSolver.expansionOrder();
                            }
                            //g0->shift(config.U/2.0);

                            //IOhelper::plot(*g0, config.beta, "Weiss Function" + std::to_string(dmftIt));

                            ioh.setIteration(dmftIt);
                            if(config.local.rank() == 0)
                            {
                                LOG(INFO) << "finished sampling. average expansion order: " << iSolver.avgN();
                                LOG(INFO) << "measuring impurity Greens function";
                            }
                            iSolver.computeImpGF();
                            //gImp->setParaMagnetic();
                            //IOhelper::plot(*gImp, config.beta, "after measure gIp Function " + std::to_string(dmftIt));
                            if(config.local.rank() == 0)
                            {
                                LOG(INFO) << "forcing paramagnetic solution";
                                LOG(INFO) << "Writing results";
                            }
                            //
                            //TODO: even for cthyb?
                            //g0*.setParaMagnetic();

                            ioh.writeToFile();

                        }

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
        };
} // end namespace DMFT
#endif
