#ifndef DMFT_BETHE_HPP_
#define DMFT_BETHE_HPP_
#include "GreensFct.hpp"
#include "WeakCoupling.hpp"
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
                DMFT_BetheLattice(const std::string& outDir, const Config& config, RealT mixing, ImpSolver &solver, GreensFct &G0, GreensFct &GImp, const RealT D):
                    config(config), mixing(mixing), iSolver(solver), g0(G0), g0Info("G0"), gImp(GImp), gImpInfo("GImp"),sImpGF(config.beta),sImpInfo("sImp"), ioh(outDir, config), D(D)
            {
                std::string tmp("G0_Guess");
                ioh.writeToFile(g0,tmp);
                ioh.addGF(g0,   g0Info);
                ioh.addGF(gImp, gImpInfo);
                ioh.addGF(sImpGF,sImpInfo);
            };

                void update(const unsigned int iterations, const unsigned long long updates){
                    for(unsigned int dmftIt = 1; dmftIt < iterations+1; dmftIt++){
                        LOG(INFO) << "Computing new Weiss Green's function";
                        //TODO: vectorize
                        for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                            for(int s=0;s<_CONFIG_spins;s++){
                                // set SC condition, enforcing Paramagnetic solution
                                int sign = (2*(mFreq(n,config.beta)>0))-1;
                                ComplexT tmp = ComplexT(config.mu,mFreq(n,config.beta)) - (D/2.0)*(D/2.0)*gImp.getByMFreq(n,s);
                                g0.setByMFreq(n,s,1.0/tmp );
                            }
                        }
                        g0.transformMtoT();
                        //TODO: this should be part of WeakCoupling
                        g0.shift(config.U/2.0);
                        g0.setParaMagnetic();
                        if(dmftIt == 1)
                        {
                            ioh.setIteration(0);
                            ioh.writeToFile();
                        }
                        
                        //only for WeakCoupling., this is now in update itself
                        //g0.shift(config.U/2.0);

                        MatG sImp = (g0.getMGF().cwiseInverse() - gImp.getMGF().cwiseInverse());
                        sImpGF.setByMFreq(sImp);
                        sImpGF.transformMtoT();

                        //TODO this breaks for large N
                        for(long unsigned int i=0; i <= 20; i++){
                            iSolver.update(updates/20.0);
                            LOG(INFO) << "U=" << config.U << ", beta=" << config.beta <<  " at iteration: " << dmftIt << "(" << (5*i) << "%) order: " << iSolver.expansionOrder();
                        }
                        //g0.shift(config.U/2.0);

                        //IOhelper::plot(g0, config.beta, "Weiss Function" + std::to_string(dmftIt));

                        LOG(INFO) << "finished sampling";
                        LOG(INFO) << "measuring impurity Greens function";
                        iSolver.computeImpGF();

                        LOG(INFO) << "forcing paramagnetic solution";
                        LOG(INFO) << "Writing results";
                        ioh.setIteration(dmftIt);
                        ioh.writeToFile();

                    }
                }


            private:
                // general settings
                const Config&	config;
                IOhelper	    ioh;

                // lattice specific
                // TODO separate class
                const RealT D;




                RealT mixing;

                ImpSolver& iSolver;
                LogInfos	g0Info;
                LogInfos	gImpInfo;
                LogInfos	sImpInfo;
                GreensFct&	g0;
                GreensFct&	gImp;
                GreensFct  sImpGF;
        };
} // end namespace DMFT
#endif
