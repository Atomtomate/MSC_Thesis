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
            FFT fft(config.beta);
            //PhysRevB.72.035122
            //MatG g0(_CONFIG_maxMatsFreq, 2);
            //ImTG g0_it(_CONFIG_maxTBins , 2);
            for(int n=-_CONFIG_maxMatsFreq/2;n<_CONFIG_maxMatsFreq/2;n++){
                const RealT mf = mFreqS(n + g.isSymmetric()*_CONFIG_maxMatsFreq/2, config.beta);
                const RealT signMF = 2*(mf>=0)-1;
                ComplexT tmp(0.0, 2.0*mf-2.0*signMF*std::sqrt(mf*mf + D*D));
                tmp /= (D*D);
                //g0(n+_CONFIG_maxMatsFreq/2, UP) = tmp;
                //g0(n+_CONFIG_maxMatsFreq/2, DOWN) = tmp;
                if(n < _CONFIG_maxMatsFreq)
                {
                    g.setByMFreq(n + g.isSymmetric()*_CONFIG_maxMatsFreq/2,UP,tmp);
                    g.setByMFreq(n + g.isSymmetric()*_CONFIG_maxMatsFreq/2,DOWN,tmp);
                }
            }
            //std::array<RealT,2> z = {.0, .0};
            //std::vector<std::array<RealT,2> > tail_coeff = {{1.0,1.0} , z, z, z};
            //fft.transformMtoT(g0, g0_it, tail_coeff, g.isSymmetric());
            //g.setByT(g0_it);
            g.fitTail();
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

        void setHybFromG0(GreensFct &g0, GreensFct &hyb, const RealT D, const Config &config)
        {
            g0.fitTail();
            const int n_min = _CONFIG_maxMatsFreq/(2*hyb.isSymmetric()) > 100 ? 100/(2*hyb.isSymmetric()) : _CONFIG_maxMatsFreq/2;
            for(int n= -_CONFIG_maxMatsFreq/2; n< _CONFIG_maxMatsFreq/2; n++)
            {
                int n_sym = n + hyb.isSymmetric()*_CONFIG_maxMatsFreq/2;
                const RealT mf = mFreqS(n_sym, config.beta);
                VLOG(6) << n << " <=> " << mf << " : " << g0.getByMFreq(n_sym,0) << " [-1]: " << ComplexT(0., mf) - 1./g0.getByMFreq(n_sym, 0);
                for(int s=0; s<_CONFIG_spins; s++)
                {
                    //if(std::abs(n) < n_min)
                    //{
                    hyb.setByMFreq(n_sym, s, ComplexT(0.,mf)  - 1./g0.getByMFreq(n_sym, s));
                    //}
                    //else
                    //{
                    //    hyb.setByMFreq(mf, s, D*D*0.25/ComplexT(0.,mf));
                    //}
                }
            }
            hyb.markMSet();
            hyb.transformMtoT();
            hyb.shiftZero();
        }

        void setSCG0(GreensFct &g0, const Config &config)
        {
            LOG(WARNING) << "simple cubic guess has numerical problems in this version";
            for(int n=0;n<_CONFIG_maxMatsFreq;n++)
            {
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

        void _run_hysteresis(RealT beta, boost::mpi::communicator world, boost::mpi::communicator local, const bool isGenerator)
        {
            const int D = 1.0;          // half bandwidth
            const RealT t	= D/2.0;
            const int burnin = 30000;
            const RealT zeroShift = 0.505;
            DMFT::GreensFct g0(beta, true, true);				// construct new Weiss green's function
            DMFT::GreensFct gImp(beta, true, true);
            DMFT::GreensFct gLoc(beta);

            const std::string solverType = "CT-INT";
            for(double U = 0.2; U < 8.2; U += 0.2)
            {   
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << "hysteresis_Uincr_b" << beta << "U" << U;
                std::string descr = stream.str();
                Config conf(beta, U/2.0, U, D, _CONFIG_maxMatsFreq, _CONFIG_maxTBins, world,local, isGenerator, "CT-INT"); 
                LOG(INFO) << "initializing: " << descr;
                setBetheSemiCirc(g0, D, conf);
                setBetheSemiCirc(gImp, D, conf);
                DMFT::WeakCoupling impSolver(&g0, &gImp, &conf, zeroShift, burnin);
                DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, conf, 0, impSolver, &g0, &gImp, D);

                LOG(INFO) << "80k MC samples";
                dmftSolver.setDir(descr + "_1M");
                dmftSolver.solve(10,80000);
                LOG(INFO) << "500k MC samples";
                dmftSolver.solve(10,500000);
                LOG(INFO) << "1M MC samples";
                dmftSolver.solve(10,1000000);
                LOG(INFO) << "5M MC samples";
                dmftSolver.setDir(descr + "_50M");
                dmftSolver.solve(5,5000000);
                LOG(INFO) << "10M MC samples";
                dmftSolver.solve(5,10000000);
                LOG(INFO) << "50M MC samples";
                dmftSolver.solve(5,50000000);
                //dmftSolver.setDir(descr + "_100M");
                //LOG(INFO) << "100M MC samples";
                //dmftSolver.solve(3,100000000);
            }

            LOG(INFO) << "Converged for U = 0.2 to U = 8.0 and beta = " << beta << ".\n Computing with U = 8.0 as initial guess now.";
            for(double U = 0.2; U < 8.2; U += 0.2)
            {   
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << "hysteresis_Udecr_b" << beta << "U" << U;
                std::string descr = stream.str();
                Config conf(beta, U/2.0, U, D, _CONFIG_maxMatsFreq, _CONFIG_maxTBins, world,local,isGenerator, solverType); 
                LOG(INFO) << "initializing: " << descr;
                DMFT::GreensFct g0_2(beta, true, true);				// construct new Weiss green's function
                DMFT::GreensFct gImp_2(beta, true, true);
                DMFT::GreensFct gLoc_2(beta);
                g0_2.setByT(gImp.getItGF());
                g0_2.setByMFreq(gImp.getMGF());
                gImp_2.setByT(gImp.getItGF());
                gImp_2.setByMFreq(gImp.getMGF());

                DMFT::WeakCoupling impSolver(&g0_2, &gImp_2, &conf, zeroShift, burnin);
                DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, conf, 0, impSolver, &g0_2, &gImp_2, D);

                LOG(INFO) << "80k MC samples";
                dmftSolver.setDir(descr + "_1M");
                dmftSolver.solve(10,80000);
                LOG(INFO) << "500k MC samples";
                dmftSolver.solve(10,500000);
                LOG(INFO) << "1M MC samples";
                dmftSolver.solve(10,1000000);
                LOG(INFO) << "5M MC samples";
                dmftSolver.setDir(descr + "_50M");
                dmftSolver.solve(5,5000000);
                LOG(INFO) << "10M MC samples";
                dmftSolver.solve(5,10000000);
                LOG(INFO) << "50M MC samples";
                dmftSolver.solve(5,50000000);

            }
        }

        void _test_hysteresis(boost::mpi::communicator world, boost::mpi::communicator local, const bool isGenerator)
        {

            _run_hysteresis(10.0, world, local, isGenerator);
            /*std::thread t1(_run_hysteresis,U0);
              std::thread t2(_run_hysteresis,U5);
              t1.join();
              t2.join();*/
        }

        // =========== 	generate input		============
        // -- atomic limit, H = t \sum_nn i->j - mu \sum_sites,spin n_i,s + U sum_sites n_i,d n_i,u
        // _test_SOH -- non interacting lmit: U=0
        // _test_SOH -- U>0

        int _test_PT(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator, bool use_bethe, double mixing)
        {

            std::string solverType = "CT-INT";
            if(isGenerator)
            {
                DMFT::Config config(64.0, 1.3, 2.6, 1.0, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType);
                if(config.world.rank() == 0){
                    LOG(INFO) << "Testing with Bethe lattice guess, sc energy density, U = " << config.U;
                    LOG(INFO) << "Setting D = 1, t = D/2.0, a = 1";
                }
                const int D = 1.0;          // half bandwidth
                const RealT a	= 1.0;
                const RealT t	= D/2.0;
                const int burnin = 30000;
                const RealT zeroShift = 0.5045;//0.5016;

                std::string descr = "BetheLatticePT";

                DMFT::GreensFct g0  (config.beta, true, true);							// construct new Weiss green's function
                DMFT::GreensFct gImp(config.beta, true, true);
                DMFT::GreensFct gLoc(config.beta);

                if(config.world.rank() == 0){
                    LOG(INFO) << "computing analytical solution for G_wn for U=0";
                    //std::cout << "G0 guess: 0 for bethe, 1 for SC: ";
                    //std::cin >> guess;
                    ////if(!use_bethe){
                    //LOG(INFO) << "using simple cubic guess";
                    //setSCG0(g0,config);
                    //}else{
                    LOG(INFO) << "using bethe lattice guess";
                }

                //TODO: write hilbert transform function
                setBetheSemiCirc(gImp, D, config);
                setBetheSemiCirc(g0, D, config);
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

                LOG(INFO) << "initializing rank " << config.world.rank() << ". isGenerator == " << config.isGenerator ;
                DMFT::WeakCoupling impSolver(&g0, &gImp, &config, zeroShift, burnin);
                if(use_bethe)
                {
                    DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config, mixing, impSolver, &g0, &gImp, D);
                    //TODO: delete old or keep it
                    LOG(INFO) << "1 Million MC steps";
                    dmftSolver.setDir("BetheLatticePT_it1_80k");
                    dmftSolver.solve(10,80000);

                    LOG(INFO) << "5 Million MC steps";
                    dmftSolver.setDir("BetheLatticePT_it2_700k");
                    dmftSolver.solve(10,700000);

                    LOG(INFO) << "10 Million MC steps";
                    dmftSolver.setDir("BetheLatticePT_it3_10m");
                    dmftSolver.solve(10,1000000);

                    LOG(INFO) << "50 Million MC steps";
                    dmftSolver.setDir("BetheLatticePT_it4_50m");
                    dmftSolver.solve(5,5000000);

                    LOG(INFO) << "100 Million MC steps";
                    dmftSolver.setDir("BetheLatticePT_it5_100m");
                    dmftSolver.solve(5,100000000);

                    LOG(INFO) << "500 Million MC steps";
                    dmftSolver.setDir("BetheLatticePT_it6_500m");
                    dmftSolver.solve(2,1000000000);

                    LOG(INFO) << "800 Million MC steps";
                    dmftSolver.setDir("BetheLatticePT_it7_800m");
                    dmftSolver.solve(2,1000000000);
                }
                else
                {
                    LOG(WARNING) << "simple cubic lattice DMFT currently disabled!";
                    //REMARK: paper for type deduction from constructor has been accepted for C++17
                    //DMFTSolver<WeakCoupling> dmftSolver(descr, config, mixing, impSolver, g0, gImp, gLoc);
                    //dmftSolver.update(30,1000000);
                }
                //Bethe solution: sImp = t_t1*t_t1*gImp.getByMFreq(n,s);

            }
            return 0;
        }

        int _test_SOH(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
        {
            const int D = 1.0;          // half bandwidth
            const RealT a	= 1.0;
            const RealT t	= D/2.0;
            const int burnin = 10000;
            const RealT zeroShift = 0.51;

            const RealT beta    = 10;

            const std::string solverType = "CT-INT";
            if(isGenerator)
            {
                for(int U_l : {1,2,3,4,5})
                {
                    std::string descr = "SBHubbardTest/U" + std::to_string(static_cast<int>(U_l));
                    const RealT U       = U_l;
                    const RealT mu      = U/2.0;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType);
                    DMFT::GreensFct g0  (config.beta, true, true);					// construct new Weiss green's function
                    DMFT::GreensFct gImp(config.beta, true, true);
                    DMFT::GreensFct gLoc(config.beta);
                    setBetheSemiCirc(g0, D, config);
                    setBetheSemiCirc(gImp, D, config);

                    LOG(INFO) << "initializing rank " << config.world.rank() << ". isGenerator ==" << config.isGenerator ;
                    DMFT::WeakCoupling impSolver(&g0, &gImp, &config, zeroShift, burnin);
                    DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config, 0.0, impSolver, &g0, &gImp, D);
                    dmftSolver.solve(5*U_l, 100000, true);

                    //(config.local.barrier)();                                           // wait to finish
                }
            }
        }

DMFT::ComplexT fit_sym_tail(int n, int i, const RealT beta)
{
    RealT wn = DMFT::mFreqS(n, beta);
    ComplexT iwn = ComplexT(0., wn);
    //if(i == -1) return iwn;
    //if(i == 1) return 1.0/iwn;
    //if(i == 3) return -1.0/(wn*wn*iwn);
    //if(i == 5) return 1./(wn*wn*wn*wn*iwn);
    if(i%2 == 1 && i > 0) return std::pow(iwn, -i);
    return ComplexT(0.,0.);
}
DMFT::ComplexT basic_tail(int n, int i, const RealT beta)
{
    if(i == 1) return 1./ComplexT(0., DMFT::mFreqS(n, beta));
    return ComplexT(0.,0.);
}
DMFT::ComplexT tmp(DMFT::ComplexT x, int i)
{
    return 1.0/std::pow(x,i);
}

        int _test_hyb(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
        {
            const int D = 2.0;          // half bandwidth
            const RealT a	= 1.0;
            const RealT t	= D/2.0;
            const int burnin = 30000;

            const RealT beta    = 64;
            const int U_l       = 10;
            GFTail tail;
            //hybTail.fitFct = &fit_sym_tail;
            tail.fitFct = &fit_sym_tail;
            tail.nC = 6;
            tail.first = 20;
            tail.last = 200 ;


            DMFT::GreensFct g0(beta, true, true, tail);					// construct new Weiss green's function
            DMFT::GreensFct* hyb = new  DMFT::GreensFct(beta, true, true, tail);					// construct new hybridization function
            DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
            DMFT::GreensFct gLoc(beta);
            const std::string solverType = "CT-HYB";
            //for(double U_l = 0; U_l < 7; U_l += 2)
            //{   
                if(isGenerator)
                {
                    std::string descr = "SBHubbardHyb_U" + std::to_string(static_cast<int>(U_l));
                    const RealT U       = U_l;
                    const RealT mu      = U/2.0;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType);
                    setBetheSemiCirc(g0, D, config);
                    setBetheSemiCirc(*gImp, D, config);
                    LOG(INFO) << "Setting hybridization function from Weiss function guess.";
                    setHybFromG0(g0, *hyb, D, config);
                    //IOhelper::plot(*hyb, config.beta, "Hybridization Function");
                    //IOhelper::plot(g0, config.beta, "Weiss Function");
                    IOhelper::plot(*gImp, config.beta, "Impurity Function");

                    LOG(INFO) << "initializing rank " << config.world.rank() << ". isGenerator ==" << config.isGenerator;
                    DMFT::StrongCoupling impSolver(hyb, gImp, &config, burnin);
                    DMFT::DMFT_BetheLattice<StrongCoupling> dmftSolver(descr, config, 0.0, impSolver, hyb, gImp, D, true);
                    dmftSolver.solve(25,8000000, true);

                    //(config.local.barrier)();                                           // wait to finish
                //
                }
            //}
            delete(gImp);
            delete(hyb);
        }

    void _test_average_PO(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
            const int D = 1.0;          // half bandwidth
            const RealT a	= 1.0;
            const RealT t	= D/2.0;
            const int burnin = 20000;
            const RealT beta  = 128;
            DMFT::GFTail tail;
            const RealT zeroShift = 0.504;
            tail.fitFct = &fit_sym_tail;
            std::string descrINT = "expOrder_CTINT";
            const std::string solverTypeI = "CT-INT";
            std::string descrHYB = "expOrder_CTHYB";
            const std::string solverTypeH = "CT-HYB";
            for(RealT U : {0.0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.0, 5.5, 6.0})//{})
            {
                DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                DMFT::GreensFct* hyb = new  DMFT::GreensFct(beta, true, true, tail);					// construct new hybridization function
                DMFT::GreensFct gLoc(beta);
                RealT mu    = U/2.;
                DMFT::Config configHYB(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverTypeH, descrHYB);
                DMFT::Config configINT(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverTypeI, descrINT);

                LOG(INFO) << "setting up initial guess for CT-INT at U = " << U;
                setBetheSemiCirc(*g0, D, configINT);
                setBetheSemiCirc(*gImp, D, configINT);
                LOG(INFO) << "initializing rank " << configINT.world.rank() << ". isGenerator ==" << configINT.isGenerator ;
                DMFT::WeakCoupling impSolverINT(g0, gImp, &configINT, zeroShift, burnin);
                DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolverINT(descrINT, configINT, 0.0, impSolverINT, g0, gImp, D);
                dmftSolverINT.solve(1, 8000000, true);


                LOG(INFO) << "setting up initial guess for CT-HYB at U = " << U;
                setBetheSemiCirc(*g0, D, configHYB);
                setBetheSemiCirc(*gImp, D, configHYB);
                setHybFromG0(*g0, *hyb, D, configHYB);
                LOG(INFO) << "initializing rank " << configHYB.world.rank() << ". isGenerator ==" << configHYB.isGenerator ;
                DMFT::StrongCoupling impSolverHYB(hyb, gImp, &configHYB, burnin);
                DMFT::DMFT_BetheLattice<StrongCoupling> dmftSolverHYB(descrHYB, configHYB, 0.0, impSolverHYB, hyb, gImp, D, true);
                dmftSolverHYB.solve(1, 8000000, true);
   
                //exit(0);
                delete(g0);
                delete(hyb);
                delete(gImp);
            }
    }

    void _test_IPT(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
            const int D = 1.0;          // half bandwidth
            const RealT t	= D/2.0;

            const RealT beta  = 40;
            DMFT::GFTail tail;
            tail.fitFct = &fit_sym_tail;
            std::string descr = "IPT_Bethe_PT";
            const std::string solverType = "IPT";
            for(RealT U : {0., 1., 2., 2.5, 3., 4.})
            {
                DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                RealT mu    = U/2.;
                DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType);
                setBetheSemiCirc(*gImp, D, config);
                auto ipt_solver = DMFT::IPT(descr, g0, gImp, config, D);
                LOG(INFO) << "Solve impurity problem using IPT for U = " + std::to_string(U);
                ipt_solver.solve(20, 100, false);
                delete(gImp);
                delete(g0);
            }
    }
    void _IPT_Z(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
            const int D = 1;          // half bandwidth
            const RealT t	= D/2.0;
            const std::string solverType = "IPT";
            DMFT::GFTail tail;
            tail.fitFct = &fit_sym_tail;
            RealT U_inc = 0.01;
            std::vector<DMFT::GreensFct> g0_bak;
            std::vector<DMFT::GreensFct> gImp_bak;
            std::array<RealT,3> betal = {133., 200., 400};
            RealT U = 0.0;
            for(int i = 0; i < betal.size(); i++)//
            {
                RealT beta = betal[i];
                g0_bak.emplace_back(beta, true, true, tail);
                gImp_bak.emplace_back(beta, true, true, tail);;
                U = 0.;
                std::string descr = "IPT_Bethe_Z_from_0__";
                //for(RealT U : {0.5, 1.0, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3., 4.0})
                while(U < 4.0) //
                {
                    DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                    RealT mu    = U/2.;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
                    setBetheSemiCirc(*gImp, D, config);
                    auto ipt_solver = DMFT::IPT(descr, g0, gImp, config, D);
                    LOG(INFO) << "Solving impurity problem using IPT for beta = " +std::to_string(beta) + ", U = " + std::to_string(U);
                    ipt_solver.solve(50, 100, false, true);
                    U += U_inc;
                    if(U>=4.0 && beta == 400.)
                    {
                        g0_bak[i] = (*g0);
                        gImp_bak[i] = (*gImp);
                    }
                    delete(gImp);
                    delete(g0);
                }
            }
            std::string descr = "IPT_Bethe_Z_to_0";
            LOG(INFO) << "Using converged solution as initial guess";
            for(int i = 0; i < betal.size(); i++)//
            {
                RealT beta = betal[i];
                while(U >= 0.)
                {
                    DMFT::GreensFct g0(g0_bak[i]);
                    DMFT::GreensFct gImp(gImp_bak[i]);
                    RealT mu    = U/2.;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
                    auto ipt_solver = DMFT::IPT(descr, &g0, &gImp, config, D);
                    LOG(INFO) << "Solving impurity problem using IPT for beta = " +std::to_string(beta) + ", U = " + std::to_string(U);
                    ipt_solver.solve(50, 100, false, true);
                    U -= U_inc;
                }
            }

    }

    void _IPT_PD(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
            const int D = 1;          // half bandwidth
            const RealT t	= D/2.0;
            const std::string solverType = "IPT";

            for(RealT beta : {5., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 30., 34., 38., 42., 44., 46., 50., 55., 60., 65., 70., 80., 90., 100., 133., 200., 400.})//
            {
                DMFT::GFTail tail;
                tail.fitFct = &fit_sym_tail;
                std::string descr = "IPT_Bethe_PD";
                //for(RealT U : {0.5, 1.0, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3., 4.0})
                for(RealT U : {2.0, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.4, 4.6, 4.8, 5.0}) //
                {
                    DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                    RealT mu    = U/2.;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
                    setBetheSemiCirc(*gImp, D, config);
                    auto ipt_solver = DMFT::IPT(descr, g0, gImp, config, D);
                    LOG(INFO) << "Solving impurity problem using IPT for beta = " +std::to_string(beta) + ", U = " + std::to_string(U);
                    ipt_solver.solve(50, 100, false, true);
                    delete(gImp);
                    delete(g0);
                }
            }
    }
    }	//end namespace examples
}	//end namespace DMFT
