#include "singleOrbHubbard.hpp"

namespace DMFT
{
    namespace examples
    {

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

    void _test_FFT(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
        const RealT D = 1.0;          // half bandwidth
        const RealT beta  = 40;
        RealT U = 1.0;
        RealT mu    = U/2.;
        DMFT::GFTail tail;
        tail.fitFct = &fit_sym_tail;
        std::string descr = "FFT_test";
        const std::string solverType = "None";
        DMFT::GreensFct g0(beta, true, true, tail);
        DMFT::GreensFct gImp(beta, true, true, tail);
        DMFT::LogInfos gImpInfo("GImp", true, true, true);
        DMFT::LogInfos g0Info("G0", true, true, true);
        DMFT::LogInfos seInfo("Sigma", true, true, true);
        DMFT::GreensFct selfE(beta, true, true);
        DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType);
        DMFT::IOhelper ioh(descr, config);
        ioh.addGF(g0,   g0Info);
        ioh.addGF(gImp, gImpInfo);
        ioh.addGF(selfE, seInfo);
        setBetheSemiCirc(gImp, D, config);
        ioh.setIteration(0);
        for(int f = 0; f < _CONFIG_spins; f++){
            for(int n = 0; n < _CONFIG_maxMatsFreq; n++) 
            {
                int n_g0 = n - (!g0.isSymmetric())*_CONFIG_maxMatsFreq/2;
                auto iwn = ComplexT(0.,mFreqS(n_g0, config.beta));
                g0.setByMFreq(n_g0, f, 1./(iwn - D*D*gImp.getByMFreq(n_g0, f)/4.0));
            }
        }
        g0.transformMtoT();
        for(int f = 0; f < _CONFIG_spins; f++){
        for(unsigned int it = 0; it < _CONFIG_maxTBins; it++)
        {
            RealT t = config.beta*(it)/(_CONFIG_maxTBins);
            selfE.setByT(t, f, config.U*config.U*g0(t,f)*g0(t,f)*g0(t,f));
        }
        }
        selfE.transformTtoM();
        const unsigned int runs = 200;
        LOG(INFO) << "FFT test prepared, running " << runs << " iterations of Matsubara <-> iTime transformations";
        ioh.writeToFile();
        for(int i=0; i < runs; i++)
        {
            g0.transformMtoT();
            g0.transformTtoM();
            gImp.transformMtoT();
            gImp.transformTtoM();
            selfE.transformMtoT();
            selfE.transformTtoM();
            ioh.setIteration(i+1);
            ioh.writeToFile();
        }
        
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
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, conf, impSolver, &g0, &gImp, D);

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
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, conf, impSolver, &g0_2, &gImp_2, D);

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
            const int D = 1.0;          // half bandwidth
            const RealT t	= D/2.0;
            const int burnin = 30000;
            const RealT zeroShift = 0.5045;//0.5016;
            const RealT beta = 64.;
            const RealT U = 2.6;
            const RealT mu = U/2.;

            DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType);
            std::string descr = "BetheLatticePT";

            if(config.world.rank() == 0){
                LOG(INFO) << "Testing with Bethe lattice guess, sc energy density, U = " << config.U;
                LOG(INFO) << "Setting D = 1, t = D/2.0";
            }
            DMFT::GreensFct g0  (config.beta, true, true);
            DMFT::GreensFct gImp(config.beta, true, true);
            DMFT::GreensFct gLoc(config.beta);

            if(config.world.rank() == 0){
                LOG(INFO) << "computing analytical solution for G_wn for U=0";
                LOG(INFO) << "using bethe lattice guess";
            }

            setBetheSemiCirc(gImp, D, config);
            setBetheSemiCirc(g0, D, config);


            LOG(INFO) << "initializing rank " << config.world.rank() << ". isGenerator == " << config.isGenerator ;
            DMFT::WeakCoupling impSolver(&g0, &gImp, &config, zeroShift, burnin);
            if(use_bethe)
            {
                DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config,impSolver, &g0, &gImp, D);
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
            }
        }
        return 0;
    }

    int _test_convergence(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
        const int D = 1.0;          // half bandwidth
        const RealT t	= D/2.0;
        const int burnin = 20000;
        const RealT zeroShift = 0.01;
        const RealT mixing = 0.7;
        GFTail tail;
        tail.fitFct = &fit_sym_tail;
        std::string descr = "test_convergence"; 
        const std::string solverType = "CT-INT";
        RealT beta = 20;
        RealT U = 2.0;
        RealT mu = U/2.;
        DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
        DMFT::GreensFct* g0_f = new DMFT::GreensFct(beta, true, true, tail);
        DMFT::GreensFct* gImp_f = new DMFT::GreensFct(beta, true, true, tail);
        setBetheSemiCirc(*g0_f, D, config);
        setBetheSemiCirc(*gImp_f, D, config);
        {
            DMFT::WeakCoupling impSolver(g0_f, gImp_f, &config, zeroShift, burnin);
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config, impSolver, g0_f, gImp_f, D, false);
            dmftSolver.solve(20, 3000000,g0_f, true, true);
            impSolver.reset();
        }
        // no mixing
        LOG(INFO) << "\n\nComparing convergence to no mixing";
        descr = "test_convergence_no_mixing"; 
        DMFT::Config config1(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
        DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
        DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
        setBetheSemiCirc(*g0, D, config1);
        setBetheSemiCirc(*gImp, D, config1);
        {
            DMFT::WeakCoupling impSolver(g0, gImp, &config1, zeroShift, burnin);
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config1, impSolver, g0, gImp, D, false, DMFT::Mixing::nothing, 0.0);
            dmftSolver.solve(15, 3000000,g0_f, true, true);
            impSolver.reset();
        }
        //mixing
        LOG(INFO) << "\n\nComparing convergence to mixing = 0.5";
        descr = "test_convergence_mixing"; 
        DMFT::Config config2(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
        g0 = new DMFT::GreensFct(beta, true, true, tail);
        gImp = new DMFT::GreensFct(beta, true, true, tail);
        setBetheSemiCirc(*g0, D, config2);
        setBetheSemiCirc(*gImp, D, config2);
        {
            DMFT::WeakCoupling impSolver(g0, gImp, &config2, zeroShift, burnin);
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config2, impSolver, g0, gImp, D, false, DMFT::Mixing::mixing, 0.5);
            dmftSolver.solve(15, 3000000,g0_f, true, true);
            impSolver.reset();
        }

        // bad broyden
        LOG(INFO) << "\n\nComparing convergence to bad broyden";
        descr = "test_convergence_bad_broyden"; 
        DMFT::Config config3(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
        g0 = new DMFT::GreensFct(beta, true, true, tail);
        gImp = new DMFT::GreensFct(beta, true, true, tail);
        setBetheSemiCirc(*g0, D, config3);
        setBetheSemiCirc(*gImp, D, config3);
        {
            DMFT::WeakCoupling impSolver(g0, gImp, &config3, zeroShift, burnin);
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config3, impSolver, g0, gImp, D, false, DMFT::Mixing::bad_broyden, 0.5);
            dmftSolver.solve(15, 3000000,g0_f, true, true);
            impSolver.reset();
        }

        // good broyden
        LOG(INFO) << "\n\nComparing convergence to good broyden";
        descr = "test_convergence_good_groyden"; 
        DMFT::Config config4(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
        g0 = new DMFT::GreensFct(beta, true, true, tail);
        gImp = new DMFT::GreensFct(beta, true, true, tail);
        setBetheSemiCirc(*g0, D, config4);
        setBetheSemiCirc(*gImp, D, config4);
        {
            DMFT::WeakCoupling impSolver(g0, gImp, &config4, zeroShift, burnin);
            DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr, config4, impSolver, g0, gImp, D, false, DMFT::Mixing::good_broyden, 0.5);
            dmftSolver.solve(15, 3000000,g0_f, true, true);
            impSolver.reset();
        }
    }


    int _test_SOH(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
        const int D = 1.0;          // half bandwidth
        const RealT t	= D/2.0;
        const int burnin = 50000;
        const RealT zeroShift = 0.01;
        const RealT mixing = 0.2;
        //const RealT beta    = 20;
        constexpr bool readf = false;
        constexpr bool run_ins = false;

        GFTail tail;
        tail.fitFct = &fit_sym_tail;

        std::string descr_IG_ME = "50M_highmf_every100_m20_CTINT_PD_IG_M"; //"35M_b50_80mf_every120_gb_CTINT_PD_IG_M"; 
        std::string descr_IG_IN = "50M_highmf_every100_m20_CTINT_PD_IG_I"; //"35M_b50_80mf_every120_gb_CTINT_PD_IG_I"; 
        std::string descrTMP = "tmp"; 
        const long int mc_steps = 30000000;
        const std::string solverType = "CT-INT";
        //std::vector Ulist = {4.0, 2.8, 2.2 ,2.4, 2.6}; //,  3.6,2.4, //2.5,2.3,
        //std::vector Ulist = {4.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8};
        //std::vector Ulist = {4.0, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.0};

        std::vector Ulist = {2.38, 2.37}; // for exp order
        if(run_ins){
            Ulist = {3.2, 2.390, 2.388, 2.392, 2.394, 2.396, 2.398, 2.402, 2.404}; //3.2, 2.38, 2.37, 2.39, 2.45, 2.40 for decreasing U
        }else{
            Ulist = {2.6, 2.5, 2.46, 2.47, 2.48, 2.49}; // beta == 80
            //Ulist = {2.402, 2.404, 2.406, 2.408, 2.412, 2.414}; //  2.41, 2.42, 2.4,2.35, 2.5,  2.3 for increasing U
            //Ulist = {  2.40,2.43, 2.42, 2.415, 2.425, 2.435,2.45, 2.5}; //  2.41, 2.42, 2.4,2.35, 2.5,  2.3 for increasing U
        }
        if(isGenerator)
        {
            for(RealT beta : {70})// 80, 90}) //20, 25, 45, 50, 60, 70,})//, 
            {
                DMFT::GreensFct* g0_bak = new DMFT::GreensFct(beta, true, true, tail);
                DMFT::GreensFct* gImp_bak = new DMFT::GreensFct(beta, true, true, tail);
                for(RealT U : Ulist)
                {
                    if((run_ins && U != 3.2) || (!run_ins && U == 3.2)) continue;
                    //if(U == 3.2) continue;
                    const RealT mu      = U/2.0;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr_IG_ME);

                    DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gLoc = new DMFT::GreensFct(beta, true, true, tail);
                    setBetheSemiCirc(*g0, D, config);
                    setBetheSemiCirc(*gImp, D, config);

                    LOG(INFO) << "running beta = " << beta << ", U = " << U << ", initial guess metal. Delta = " << zeroShift << ", MC steps = " << mc_steps;
                    /*{
                        DMFT::WeakCoupling impSolver(g0, gImp, &config, zeroShift, burnin);
                        DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descrTMP, config, impSolver, g0, gImp, D, false);
                        dmftSolver.solve(20, 700000, true, false);
                        impSolver.reset();
                    }
                    {
                        DMFT::WeakCoupling impSolver(g0, gImp, &config, zeroShift, burnin);
                        DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descrTMP, config, impSolver, g0, gImp, D, false);
                        dmftSolver.solve(10, 1000000, true, false);
                        impSolver.reset();
                    }*/
                    {
                        DMFT::WeakCoupling impSolver(g0, gImp, &config, zeroShift, burnin);
                        DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr_IG_ME, config,impSolver, g0, gImp, D, false, DMFT::Mixing::mixing, mixing);
                        dmftSolver.solve(50, mc_steps, true, true);
                        impSolver.reset();
                    }
                    //(config.local.barrier)();
                    if(U == 3.2)
                    {
                        (*g0_bak) = (*g0);
                        (*gImp_bak) = (*gImp);
                    }
                    delete(g0);
                    delete(gImp);
                    delete(gLoc);
                }
                for(RealT U : Ulist)
                {
                    if(!run_ins) continue;
                    const RealT mu      = U/2.0;
                    DMFT::Config config2(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr_IG_IN);
                    DMFT::GreensFct* g02 =  new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gImp2 =  new DMFT::GreensFct(beta, true, true, tail);
                    if(readf)
                    {
                        std::string tmp = descr_IG_ME; //
                        std::string tmp2 = "f_GImp_b50_000000_U3_200000.out";
                        DMFT::IOhelper ioh(tmp, config2);
                        ioh.readFromFile(*gImp2, tmp2,tmp2);

                        CVectorT G0_out(_CONFIG_maxMatsFreq);
                        for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                            const int n_g0 = n + ((int)g02->isSymmetric() - 1)*_CONFIG_maxMatsFreq/2;
                            ComplexT tmp = ComplexT(0., mFreqS(n_g0, config2.beta)) - (D/2.0)*(D/2.0)*gImp2->getByMFreq(n_g0, 0);
                            G0_out(n) = 1./(tmp);
                        }
                        G0_out.real() = G0_out.real()*0.;
                        g02->g_wn.col(0) = G0_out;
                        g02->g_wn.col(1) = G0_out;
                        g02->markMSet();g02->transformMtoT();
                    }else{
                        (*g02) =    (*g0_bak);
                        (*gImp2) =    (*gImp_bak);
                    }
                    if(U == 3.2) continue;


                    LOG(INFO) << "running beta = " << beta << ", U = " << U << ", initial guess insulator";
                    {
                        DMFT::WeakCoupling impSolver2(g02, gImp2, &config2, zeroShift, burnin);
                        DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolver(descr_IG_IN, config2, impSolver2, g02, gImp2, D, false, DMFT::Mixing::mixing, mixing);
                        dmftSolver.solve(50, mc_steps, true, true);
                        impSolver2.reset();
                    }
                    //(config.local.barrier)();
                    delete(g02);
                    delete(gImp2);
                }
            }
        }
    }

    int _test_hyb(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
        const int D = 1.0;          // half bandwidth
        const RealT t	= D/2.0;
        const int burnin = 10000;

        const RealT beta    = 40;
        const int U_l       = 3;
        GFTail tail;
        tail.fitFct = &fit_sym_tail;
        //hybTail.fitFct = &fit_sym_tail;
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
                IOhelper::plot(*gImp, config.beta, "Impurity Function");

                LOG(INFO) << "initializing rank " << config.world.rank() << ". isGenerator ==" << config.isGenerator;
                DMFT::StrongCoupling impSolver(hyb, gImp, &config, burnin);
                DMFT::DMFT_BetheLattice<StrongCoupling> dmftSolver(descr, config, impSolver, hyb, gImp, D, true);
                dmftSolver.solve(5,300, true);

                //(config.local.barrier)();*                                           // wait to finish
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
                DMFT::DMFT_BetheLattice<WeakCoupling> dmftSolverINT(descrINT, configINT, impSolverINT, g0, gImp, D);
                dmftSolverINT.solve(1, 8000000, true);


                LOG(INFO) << "setting up initial guess for CT-HYB at U = " << U;
                setBetheSemiCirc(*g0, D, configHYB);
                setBetheSemiCirc(*gImp, D, configHYB);
                setHybFromG0(*g0, *hyb, D, configHYB);
                LOG(INFO) << "initializing rank " << configHYB.world.rank() << ". isGenerator ==" << configHYB.isGenerator ;
                DMFT::StrongCoupling impSolverHYB(hyb, gImp, &configHYB, burnin);
                DMFT::DMFT_BetheLattice<StrongCoupling> dmftSolverHYB(descrHYB, configHYB, impSolverHYB, hyb, gImp, D, true);
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
            RealT U = 2.0;
            for(int i = 0; i < betal.size(); i++)//
            {
                RealT beta = betal[i];
                g0_bak.emplace_back(beta, true, true, tail);
                gImp_bak.emplace_back(beta, true, true, tail);;
                U = 2.5;
                std::string descr = "IPT_Bethe_Z_from_0__";
                //for(RealT U : {0.5, 1.0, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3., 4.0})
                while(U <= 3.5) //
                {
                    DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                    RealT mu    = U/2.;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
                    setBetheSemiCirc(*gImp, D, config);
                    auto ipt_solver = DMFT::IPT(descr, g0, gImp, config, D);
                    LOG(INFO) << "Solving impurity problem using IPT for beta = " +std::to_string(beta) + ", U = " + std::to_string(U);
                    ipt_solver.solve(80, 400, false, true);
                    U += U_inc;
                    if(U>=2.5 && beta == 400.)
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
                U = 3.5;
                while(U >= 2.5)
                {
                    DMFT::GreensFct g0(g0_bak[i]);
                    DMFT::GreensFct gImp(gImp_bak[i]);
                    RealT mu    = U/2.;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
                    auto ipt_solver = DMFT::IPT(descr, &g0, &gImp, config, D);
                    LOG(INFO) << "Solving impurity problem using IPT for beta = " +std::to_string(beta) + ", U = " + std::to_string(U);
                    ipt_solver.solve(80, 400, false, true);
                    U -= U_inc;
                }
            }

    }

    void _IPT_PD(const boost::mpi::communicator local, const boost::mpi::communicator world, const bool isGenerator)
    {
            const int D = 1;          // half bandwidth
            const RealT t	= D/2.0;
            const std::string solverType = "IPT";
            DMFT::GFTail tail;
            tail.fitFct = &fit_sym_tail;
            const RealT U_min = 2.3;
            const RealT U_max = 3.5;
            const RealT U_inc = 0.01;
            std::vector<RealT> betal(120);
            RealT T_max = 0.06;
            RealT T_min = 0.005; 
            RealT dt = (T_max - T_min)/120;
            for(int i = 0; i < 120; i++)
                betal[i] = 1./(T_min + i*dt);
            RealT U = U_min;
            for(int i = 0; i < 15; i++)//
            {
                RealT beta = betal[i];
                LOG(INFO) << "starting hysteresis for beta = " << beta;
                DMFT::GreensFct* gImp_bak = new DMFT::GreensFct(beta, true, true, tail);
                U = U_min;
                //if(i < 32) U = U_max;
                std::string descr = "IPT_PD_Z_03";
                while(U <= U_max) //
                {
                    DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                    RealT mu    = U/2.;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
                    setBetheSemiCirc(*gImp, D, config);
                    setBetheSemiCirc(*g0, D, config);
                    auto ipt_solver = DMFT::IPT(descr, g0, gImp, config, D);
                    LOG(INFO) << "Solving impurity problem using IPT for beta = " +std::to_string(beta) + ", U = " + std::to_string(U);
                    ipt_solver.solve(800, 100, false, false);
                    U += U_inc;
                    //if((U>=(U_max - U_inc)))
                    //{
                    //    (*gImp_bak) = (*gImp);
                    //}else{
                        delete(gImp);
                        delete(g0);
                    //}
                }
                std::string descrB = "IPT_PD_Z_30";
                std::string descTMP = "tmp";
                DMFT::GreensFct g0a(beta, true, true, tail);
                DMFT::GreensFct gImpa(beta, true, true, tail);
                DMFT::Config config(beta, 2., 4., D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descr);
                setBetheSemiCirc(gImpa, D, config);
                setBetheSemiCirc(g0a, D, config);
                auto ipt_solver = DMFT::IPT(descTMP, &g0a, &gImpa, config, D);
                ipt_solver.solve(1000, 0, false, false);
                (*gImp_bak) = gImpa;
                U = U_max;
                while(U >= U_min)
                {
                    DMFT::GreensFct* g0 = new DMFT::GreensFct(beta, true, true, tail);
                    DMFT::GreensFct* gImp = new DMFT::GreensFct(beta, true, true, tail);
                    (*gImp) = (*gImp_bak);
                    RealT mu    = U/2.;
                    DMFT::Config config(beta, mu, U, D, DMFT::_CONFIG_maxMatsFreq, DMFT::_CONFIG_maxTBins, local, world, isGenerator, solverType, descrB);
                    auto ipt_solver = DMFT::IPT(descrB, g0, gImp, config, D);
                    LOG(INFO) << "Solving impurity problem using IPT for beta = " +std::to_string(beta) + ", U = " + std::to_string(U) << " with IG in MI phase";
                    ipt_solver.solve(800, 100, false, false);
                    U -= U_inc;
                    delete(gImp);
                    delete(g0);
                }
                delete(gImp_bak);
            }
        }
    }	//end namespace examples
}	//end namespace DMFT
