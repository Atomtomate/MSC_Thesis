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
                DMFT_BetheLattice(std::string& outDir, const Config& config, RealT mixing, ImpSolver &solver, GreensFct &G0, GreensFct &GImp, const RealT D):
                    config(config), mixing(mixing), iSolver(solver), \
                    g0(G0), g0Info("G0"), gImp(GImp), gImpInfo("GImp"), selfE(config.beta, true, false), seLInfo("SelfE"),\
                    ioh(outDir, config), D(D), fft(config.beta)
                {
                    std::string tmp("G0_Guess");
                    ioh.writeToFile(g0,tmp);
                    ioh.addGF(g0,   g0Info);
                    ioh.addGF(gImp, gImpInfo);
                    ioh.addGF(selfE, seLInfo);
                }

                void solve(const unsigned int iterations, const unsigned long long updates, bool symmetricG0 = false)
                {
                    if(config.isGenerator)
                    {
                        //TODO: move IO from generator to accumulators
                        for(unsigned int dmftIt = 1;dmftIt < iterations+1; dmftIt++)
                        {  
                            if(config.local.rank() == 0) LOG(INFO) << "Computing new Weiss Green's function";
                            //TODO: vectorize
                            MatG tmpG0(_CONFIG_maxTBins,2);
                            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                                for(int s=0;s<_CONFIG_spins;s++){
                                    // set SC condition, enforcing Paramagnetic solution
                                    // use symmetry here
                                    ComplexT tmp = ComplexT(config.mu, mFreqS(n,config.beta)) - (D/2.0)*(D/2.0)*gImp.getByMFreq(n,s);
                                    VLOG(5) << n << ": " << ComplexT(config.mu, mFreqS(n,config.beta)) << " - " << (D/2.0)*(D/2.0)*gImp.getByMFreq(n,s) << " = " << tmp;
                                    g0.setByMFreq(n,s, 1.0/tmp );
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

//TODO: use tail
                            MatG sImp = (g0.getMGF().cwiseInverse() - gImp.getMGF().cwiseInverse());
                            ImTG sImp_it(_CONFIG_maxTBins, 2);

                            selfE.setByMFreq(sImp);
                            selfE.setParaMagnetic();

                            fft.transformMtoT(sImp,sImp_it,true); 
                            //TODO improve this naive loop
                            for(long unsigned int i=0; i <= 20; i++){
                                iSolver.update(updates/20.0);
                                LOG(INFO) << "MC Walker [" << config.world.rank() << "] at "<< " (" << (5*i) << "%) of iteration " << dmftIt << ". expansion order: " << iSolver.expansionOrder();
                            }
                            //g0.shift(config.U/2.0);

                            //IOhelper::plot(g0, config.beta, "Weiss Function" + std::to_string(dmftIt));

                            ioh.setIteration(dmftIt);
                            if(config.local.rank() == 0)
                            {
                                LOG(INFO) << "finished sampling. average expansion order: " << iSolver.avgN();
                                LOG(INFO) << "measuring impurity Greens function";
                            }
                            iSolver.computeImpGF();
                            if(config.local.rank() == 0)
                            {
                                LOG(INFO) << "forcing paramagnetic solution";
                                LOG(INFO) << "Writing results";
                            }
                            gImp.setParaMagnetic();
                            g0.setParaMagnetic();

                            ioh.writeToFile();

                        }

                        //config.world.send(0, static_cast<int>(MPI_MSG_TAGS::FINALIZE), 1 );
                    }
                    else
                    {
                        //MCAccumulator<_CONFIG_maxSBins> *mcAcc = new MCAccumulator<_CONFIG_maxSBins>(g0,config);
                        //mcAcc->collect();
                        //TODO: compute GImp
                        //TODO: compute new g0
                        //TODO: distribute g0 to generators
                    }

                }

                void setDir(std::string d) { ioh.initDir(d); }

            private:
                // general settings
                const Config&	config;
                IOhelper	    ioh;

                FFT fft;

                // lattice specific
                // TODO separate class
                const RealT D;




                RealT mixing;

                ImpSolver& iSolver;
                LogInfos	g0Info;
                LogInfos	gImpInfo;
                LogInfos	seLInfo;
                GreensFct&	g0;
                GreensFct&	gImp;
                GreensFct       selfE;
        };
} // end namespace DMFT
#endif
