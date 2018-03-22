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

enum Mixing {nothing = 0, mixing, bad_broyden, good_broyden };

template<class ImpSolver>
class DMFT_BetheLattice
{
    //"GImp_b"+std::to_string(config.beta)+"_U" + std::to_string(config.U)   _b"+std::to_string(config.beta)+"_U" + std::to_string(config.U)
    public:
        DMFT_BetheLattice(std::string& outDir, const Config config, ImpSolver& solver, GreensFct* G0, GreensFct * GImp, const RealT D, bool useHyb = false, Mixing mixing = Mixing::nothing, const RealT start_mixing = 0.):
            config(config), mixing(mixing), start_mixing(start_mixing), iSolver(solver), useHyb(useHyb), useBroyden(useBroyden), \
            g0(G0), g0Info( useHyb ? "Hyb" : "G0" ), gImp(GImp), gImpInfo("GImp"), selfE(config.beta, true, true), seLInfo("SelfE"),\
            ioh(outDir, config), D(D)
        {
            std::string tmp( useHyb ? "Hyb_Guess" : "G0_Guess");
            ioh.writeToFile(*g0,tmp);
            ioh.addGF(*g0,   g0Info);
            ioh.addGF(*gImp, gImpInfo);
            ioh.addGF(selfE, seLInfo);
            wn_grid = MatG::Zero(_CONFIG_maxMatsFreq, _CONFIG_spins);
            for(int n = 0; n < _CONFIG_maxMatsFreq; n++)
            {
                wn_grid(n,0) = ComplexT(0., mFreqS(n,config.beta));
                wn_grid(n,1) = ComplexT(0., mFreqS(n,config.beta));
            }
        }

        void solve(const unsigned int iterations, const unsigned long long avg_updates, GreensFct* compareTo, bool symmetricG0 = true, bool writeOut = true) {}
        void solve(const unsigned int iterations, const unsigned long long avg_updates, bool symmetricG0 = true, bool writeOut = true)
        {
            if(config.isGenerator)
            {

                CMatrixT iJ = -start_mixing*Eigen::MatrixXcd::Identity(_CONFIG_maxMatsFreq, _CONFIG_maxMatsFreq);
                bool converged = false;
                RealT distance = 999.;
                CVectorT old_F = g0->g_wn.col(0);
                CVectorT oldG0_in = g0->g_wn.col(0);
                RealT conv_dist = 99;
                for(unsigned int dmftIt = 1; dmftIt < iterations+1; dmftIt++)
                {  
                    ioh.setIteration(dmftIt);
                    iSolver.reset();
                    RealT updates = avg_updates;
                    /*if(conv_dist > 1 && (updates <= avg_updates))
                    {
                        updates = avg_updates;
                    }
                    else if((conv_dist < 0.5) && conv_dist > 0.05 && (updates <= 3*avg_updates))
                    {
                        updates = 3*avg_updates;
                    }
                    else if(conv_dist < 0.005 && dmftIt > 2)
                    {
                        updates = 5*avg_updates;
                    }*/
                    for(long unsigned int i=1; i <= 5; i++){
                        iSolver.update(updates/5.0);
                        LOG(INFO) << "MC Walker [" << config.world.rank() << "] at "<< " (" << (20*i) << "%) of iteration " << dmftIt << ". expansion order: " << iSolver.expansionOrder();
                    }
                    if(config.local.rank() == 0)
                    {
                        //LOG(INFO) << "finished sampling. average expansion order: " << iSolver.expansionOrder();
                        LOG(INFO) << "measuring impurity Green's function";
                    }

                    DMFT::GreensFct* gImp_bak = new DMFT::GreensFct(config.beta, true, true);
                    (*gImp_bak) = (*gImp);
                    iSolver.computeImpGF();
                    iSolver.avgN().writeResults(ioh, "");
                    gImp->setParaMagnetic();
                    if(config.local.rank() == 0)
                    {
                        //LOG(INFO) << "forcing paramagnetic solution and calculating new input";
                    }
                    MatG sImp(_CONFIG_maxMatsFreq, _CONFIG_spins);
                    sImp = (g0->getMGF().cwiseInverse() - gImp->getMGF().cwiseInverse());
                    //sImp = (wn_grid - (D/2.)*(D/2.)*gImp->getMGF() - gImp->getMGF().cwiseInverse() );
                    selfE.setByMFreq(sImp);
                    selfE.transformMtoT();
                    // DMFT equation
                    CVectorT G0_out(_CONFIG_maxMatsFreq);
                    CVectorT G0_in = g0->g_wn.col(0);
                    for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                        const int n_g0 = n + ((int)g0->isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                        //for(int s=0;s<_CONFIG_spins;s++){ //no need for distinction in paramagnetic solution
                            ComplexT tmp = ComplexT(0., mFreqS(n_g0, config.beta)) - (D/2.0)*(D/2.0)*gImp->getByMFreq(n_g0, 0);
                            //LOG(ERROR) << tmp << "=" << ComplexT(0., mFreqS(n_g0, config.beta)) << "-" << (D/2.0)*(D/2.0)*gImp->getByMFreq(n_g0, 0);
                            G0_out(n) = 1./(tmp); // force paramagnetic solution
                        //}
                    }
                     
                    G0_out.real() = G0_out.real()*0.;
                    //Broyden's method. for normal mixing, deactivate iJ update
                    CVectorT new_F = G0_out - G0_in;
                    CVectorT deltaG0 = (G0_in - oldG0_in); //.imag();
                    if(dmftIt > 1){
                        CVectorT delta_F = (new_F - old_F);//.imag();
                        distance = delta_F.norm();
                        //good:
                        //iJ = iJ + ((deltaG0 - iJ*delta_F)/(deltaG0.transpose()*iJ*delta_F))*(deltaG0.transpose()*iJ);
                        if(mixing == 3)
                        {
                            iJ = iJ + ((deltaG0 - iJ*delta_F)/(deltaG0.conjugate().transpose()*iJ*delta_F))*(deltaG0.conjugate().transpose()*iJ);
                        }
                        if(mixing == 2)
                            iJ = iJ + ((deltaG0 - iJ*delta_F)/(distance*distance))*(delta_F.conjugate().transpose());  // for complex use conjugate().transpose()
                    }
                    if(mixing != Mixing::nothing) 
                    {
                        G0_out = G0_in - iJ*new_F;
                    }
            
                    // set paramagnetic
                    g0->g_wn.col(0) = G0_out;
                    g0->g_wn.col(1) = G0_out;
                    old_F = new_F;
                    oldG0_in = G0_in;
                    // end of broyden's method
                    g0->markMSet();
                    g0->transformMtoT();
                    conv_dist = deltaG0.head(50).norm();
                    LOG(INFO) << "distance: " << conv_dist;
                    if(conv_dist < 0.0010 && dmftIt > 3)
                    {
                        converged = true;
                        LOG(INFO) << "converged after iteration " << dmftIt;
                    }
                    if(false && !useHyb)
                    {
                        g0->shift(config.U/2.0);
                        g0->markMSet();
                        g0->transformMtoT();
                    }
                    
                    //CVectorT tmp_dif = compareTo->g_wn.col(1)-g0->g_wn.col(1);
                    //    LOG(INFO) << "Absolute convergence: " << tmp_dif.norm();

                    if(writeOut)
                    {
                        //LOG(INFO) << "Writing results";
                        ioh.writeToFile();
                    }
                    if((dmftIt > iterations-3 ) || converged)
                    {
                        statAcc(*gImp);
                        statAccSE(selfE);
                    }
                    if(converged) break;
                }
                statAcc.setGF(*gImp);
                statAccSE.setGF(selfE);
                LogInfos finalGI_info("GImp_b"+std::to_string(config.beta)+"_U" + std::to_string(config.U));
                LogInfos finalSE_info("SelfE_b" +std::to_string(config.beta)+"_U" + std::to_string(config.U));
                ioh.writeFinalToFile(*gImp, finalGI_info, false, config.U, false);
                ioh.writeFinalToFile(selfE, finalSE_info, true, config.U, false);

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
        const bool useBroyden;
        Mixing mixing;
        const RealT start_mixing;

        MatG wn_grid;
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
