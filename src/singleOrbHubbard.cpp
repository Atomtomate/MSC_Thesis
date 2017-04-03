#include "singleOrbHubbard.hpp"

namespace DMFT
{
    namespace examples
    {

        /*! Set Weiss function as the inverse of the Hilbert transform for the semi-circular DOS with half bandwidth 2t.
         *  G0^{-1} = Sigma + 1/DH(i w_n + mu - Sigma), Sigma = 0 for initial guess.
         *  DH(zeta) = 2*(zeta -  Sign(Im(zeta)) * sqrt( zeta^2 - (2t)^2 )/( (2t)^2 )
         *
         *  @param[out] g0     Weiss function to be initialized
         *  @param[in]  D      Half band-width
         *  @param[in]  p      parameters for this calculation
         */
        void setBetheSemiCirc(GreensFct &g, const RealT D, const Config &config)
        {
            //PhysRevB.72.035122
            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                const RealT mf = mFreq(n, config.beta);
                const RealT signMF = 2*(mf>=0)-1;
                ComplexT tmp(0.0, 2.0*mf-2.0*signMF*std::sqrt(mf*mf + D*D));
                g.setByMFreq(n,0,tmp/(D*D));
                g.setByMFreq(n,1,tmp/(D*D));
            }
            g.markMSet();
            g.transformMtoT();
        }

        void setBetheSemiCirc_naive(GreensFct &g, const RealT D, const Config &config)
        {
            for(int n=0;n<_CONFIG_maxMatsFreq;n++)
            {
                const RealT mf = 3.141592653*(n-static_cast<int>(_CONFIG_maxMatsFreq/2) + 1)/config.beta;
                const RealT signMF = 2*(mf>=0)-1;
                ComplexT tmp(0.0, mf - signMF*std::sqrt(mf*mf + D*D) );
                tmp = 2.0*tmp/(D*D);
                g.setByMFreq(n,0,tmp);
                g.setByMFreq(n,1,tmp);
            }
            g.markMSet();
            g.transformMtoT();
        }

        void setSCG0(GreensFct &g0, const Config &config)
        {
            LOG(WARNING) << "simple cubic guess has numerical problems in this version";
            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                const RealT mf = mFreq(n, config.beta);
                const int s = 2*(mf > 0) - 1;
                ComplexT z(config.mu,mf);
                ComplexT res(0.0, -(s*boost::math::constants::pi<RealT>())*std::exp(-std::real(z*z)));
                res = res*Faddeeva::erfc(-ComplexT(0.0,1)*static_cast<RealT>(s)*z);
                if(std::isnan(std::imag(res)) or std::isinf(std::imag(res))) res = -ComplexT(0.0,1.0/mf);
                g0.setByMFreq(n,0,ComplexT(0.0, std::imag(res)));
                g0.setByMFreq(n,1,ComplexT(0.0, std::imag(res)));
            }
            g0.markMSet();
            g0.transformMtoT();

        }


        void _run_hysteresis(RealT U){
            const int D = 1.0;          // half bandwidth
            const RealT t	= D/2.0;
            std::string descr = "hysteresis" + std::to_string(static_cast<int>(U)) + "converging";
            std::string descrPT = "hysteresis" + std::to_string(static_cast<int>(U)) + "PT";
            const double beta = 64;
            Config conf(beta, U/2.0, U, _CONFIG_maxMatsFreq); 
            //Config cU6(beta, 3.0, 6.0, _CONFIG_maxMatsFreq);
            Config confPT(beta, 1.2, 2.4, _CONFIG_maxMatsFreq);
            //Config cPT_U6(beta, 1.2, 2.4, _CONFIG_maxMatsFreq);
            LOG(INFO) << "Using Bethe Lattice (guess and SC condition), computing Gimp for U = " << U << " at beta = 64";
            LOG(INFO) << "Using both converged solution to find solutions near phase transition (U = 2.4)";
            DMFT::GreensFct g0(beta);				// construct new Weiss green's function
            DMFT::GreensFct gImp(beta);
            DMFT::GreensFct gLoc(beta);
            setBetheSemiCirc(gImp, D, conf);


            //TODO: overload operators
            for(int s=0;s<_CONFIG_spins;s++){
                for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                    ComplexT tmp = (1.0/( ComplexT(0.0 , mFreq(n,beta)) - (D/2.0)*(D/2.0)*gImp.getByMFreq(n,s) ));
                    g0.setByMFreq(n,s, tmp);
                }
            }
            g0.transformMtoT();
            LOG(INFO) << "initializing";
            DMFT::WeakCoupling impSolver(g0, gImp, conf, 0.1, 100000);
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, conf, 0, impSolver, g0, gImp, D);

            LOG(INFO) << "solving";
            dmftSolver.update(5*(static_cast<int>(U)+1),5000000);

            LOG(INFO) << "Using converged solutions to obtain U=2.4 solution";
            DMFT::WeakCoupling impSolver_2(g0, gImp, confPT, 0.1, 100000);
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolverPT(descrPT, confPT, 0, impSolver_2, g0, gImp, D);
            LOG(INFO) << "solving, starting with converged solution";
            dmftSolverPT.update(30,5000000);
        }

        void _test_hysteresis(void)
        {
            double U0 = 0.1;
            double U5 = 5.0;
            std::thread t1(_run_hysteresis,U0);
            std::thread t2(_run_hysteresis,U5);
            t1.join();
            t2.join();
        }

        // =========== 	generate input		============
        // -- atomic limit, H = t \sum_nn i->j - mu \sum_sites,spin n_i,s + U sum_sites n_i,d n_i,u
        // _test_SOH -- non interacting lmit: U=0
        // _test_SOH -- U>0

        int _test_SOH(DMFT::Config& config, bool use_bethe, double mixing){


            LOG(INFO) << "Testing with Bethe lattice guess, sc energy density, U = " << config.U;
            LOG(INFO) << "Setting D = 1, t = D/2.0, a = 1";
            const int D = 1.0;          // half bandwidth
            const RealT a	= 1.0;
            const RealT t	= D/2.0;

            std::string descr = "BetheLattice";

            DMFT::GreensFct g0  (config.beta);							// construct new Weiss green's function
            DMFT::GreensFct gImp(config.beta);
            DMFT::GreensFct gLoc(config.beta);

            LOG(INFO) << "computing analytical solution for G_wn for U=0";
            //std::cout << "G0 guess: 0 for bethe, 1 for SC: ";
            //std::cin >> guess;
            ////if(!use_bethe){
            //LOG(INFO) << "using simple cubic guess";
            //setSCG0(g0,config);
            //}else{
            LOG(INFO) << "using bethe lattice guess";
            setBetheSemiCirc(gImp, D, config);
            gImp.transformMtoT();
            //TODO: write hilbert transform function
            setBetheSemiCirc(g0, D, config);
            g0.transformMtoT();
            //IOhelper::plot(g0, config.beta, "Weiss Function");

            //}
            //TODO: overload operators
            /*for(int s=0;s<_CONFIG_spins;s++){
                for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                    ComplexT tmp = (1.0/( ComplexT(config.mu - config.U/2.0, mFreq(n,config.beta)) - (D/2.0)*(D/2.0)*g0.getByMFreq(n,s) ));
                    gImp.setByMFreq(n,s, tmp);
                }
            }
            gImp.transformMtoT();
            IOhelper::plot(gImp, config.beta, "gImp");
            */

            LOG(INFO) << "initializing";
            DMFT::WeakCoupling impSolver(g0, gImp, config, 0.1, 200000);
            if(use_bethe)
            {

                DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config, mixing, impSolver, g0, gImp, D);
                dmftSolver.update(20,500000);
            }
            else
            {
                LOG(WARNING) << "simple cubic lattice DMFT currently disabled!";
                //REMARK: paper for type deduction from constructor has been accepted for C++17
                //DMFTSolver<WeakCoupling> dmftSolver(descr, config, mixing, impSolver, g0, gImp, gLoc);
                //dmftSolver.update(30,1000000);
            }



            //Bethe solution: sImp = t_t1*t_t1*gImp.getByMFreq(n,s);

            return 0;
        }

    }	//end namespace examples
}	//end namespace DMFT
