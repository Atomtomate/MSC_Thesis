#include "WeakCoupling.hpp"
//TODO: use BOOST_FOR and BOOST_MPI
//TODO: use improved estimator (p. 27 CT-INT numerics pdf) - transpose M
namespace DMFT
{

    WeakCoupling::WeakCoupling(
            GreensFct *const g0, GreensFct * const gImp, const Config *const config, const RealT zeroShift, const unsigned int burninSteps):
        config(config), g0(g0), gImp(gImp), zeroShift(zeroShift), burninSteps(burninSteps), expOrdAcc(*config)
    {
        n 		= 0;
        steps 		= 0;
        totalSign	= 0;
        currentSign	= 1;
        r_time.split(5, 0);
        r_spin.split(5, 1);
        r_insert.split(5, 2);
        r_accept.split(5, 3);
        r_shift.split(5, 4);
        r_time.split(config->local.size() , config->local.rank());           // choose sub−stream no. rank out of size streams
        r_spin.split(config->local.size() , config->local.rank());           // choose sub−stream no. rank out of size streams
        r_insert.split(config->local.size() , config->local.rank());         // choose sub−stream no. rank out of size streams
        r_accept.split(config->local.size() , config->local.rank());         // choose sub−stream no. rank out of size streams
        r_shift.split(config->local.size() , config->local.rank());         // choose sub−stream no. rank out of size streams
    }


    WeakCoupling::~WeakCoupling()
    {
    }

    void WeakCoupling::writeExpOrder(IOhelper ioh)
    {
        expOrdAcc.writeResults(ioh);
    }

    void WeakCoupling::computeImpGF(void)
    {
        RealT binCount = itBinsUP.size();
        RealT itCount = gImp->getItGF().rows();
        RealT mfCount = gImp->getMGF().rows();
        ImTG g_it = ImTG::Zero(itCount, 2);
        MatG g_wn = MatG::Zero(mfCount, 2);
        RealT intFac = config->beta/binCount;
        for(long j = 0; j< binCount; j++)
        {
            RealT tp = j*config->beta/(binCount);
            const RealT bValUP = boost::accumulators::sum(itBinsUP[j]);
            const RealT bValDOWN = boost::accumulators::sum(itBinsDOWN[j]);
            if(bValUP)
            {
#ifdef MATSUBARA_MEASUREMENT
                for(long k=0; k<mfCount;k++)
                {
                    RealT mfreq = gImp->isSymmetric() ? mFreqS(k,config->beta) : mFreq(k, config->beta);
                    g_wn(k,UP) += std::exp(ComplexT(0.0,mfreq*tp))*bValUP;
                }
#else
                for(long i=0; i<itCount; i++)
                {
                    RealT t = i*config->beta/(itCount);
                    g_it(i,UP) += (*g0)(t-tp,UP)*bValUP;
                }
#endif
            }
            if(bValDOWN)
            {
#ifdef MATSUBARA_MEASUREMENT
                for(long k=0; k<mfCount;k++)
                {
                    RealT mfreq = gImp->isSymmetric() ? mFreqS(k,config->beta) : mFreq(k, config->beta);
                    g_wn(k,DOWN) += std::exp(ComplexT(0.0,mfreq*tp))*bValDOWN;
                }
#else
                for(long i=0; i<itCount; i++)
                {
                    RealT t = i*config->beta/(itCount);
                    g_it(i,DOWN) += (*g0)(t-tp,DOWN)*bValDOWN;
                }
#endif
            }
        }
#ifdef MATSUBARA_MEASUREMENT
        g_wn = g0->getMGF() - g0->getMGF()*g_wn/totalSign;
        gImp->setByMFreq(g_wn);
        gImp->transformMtoT();
#else
        g_it = g0->getItGF() - g_it/totalSign;
        gImp->setByT(g_it);
        gImp->transformTtoM();
#endif
    }

    void WeakCoupling::reset()
    {
        n = 0;
        totalSign = 0;
        currentSign = 1;
        steps = 0;
        MatrixT tmp;
        VectorT tmp2;
        for(int s = 0; s< _CONFIG_maxSBins; s++)
        {
            itBinsUP[s] = AccT();
            itBinsDOWN[s] = AccT();
        }
        M[UP] = tmp;
        M[DOWN] = tmp;
        confs.clear();
        //gfCache[UP] = tmp2;
        //gfCache[DOWN] = tmp2;
        expOrdAcc.reset();
    }

    //TODO:	do this properly: accumulator, overflow, kahan, SBin not *beta, mats != tbins 
    void WeakCoupling::computeImpGF_OLD(void)
    {
        (config->local.barrier)();                                           // wait for all MC runners to finish
        //LOG(INFO) << "rank: " << config->local.rank() << " in computeImpGF";
        //if(config->local.rank() == 0)
        config->world.send(0, static_cast<int>(MPI_MSG_TAGS::SAMPLING_END));     // inform accumulators, that we are done

        std::vector<RealT> mpi_gImp;
        boost::mpi::broadcast(config->world, mpi_gImp, 0);
        ImTG gImpTMP = Eigen::Map<ImTG> (&mpi_gImp[0], mpi_gImp.size()/2  , 2);
        gImp->setByT(gImpTMP); 
        std::vector<RealT> tmpIT(_CONFIG_maxTBins,0);
        /*for(int i=0;i<_CONFIG_maxTBins;i++)
          tmpIT[i] = gImp->getItGF()(i,DOWN);
          LOG(INFO) << "----------------- mcAcc ------------------: " << totalSign << "\n"
          << tmpIT << "\n"
          << "------------------> mcAcc <---------------------";
          */
        gImp->transformTtoM();
    }


    void WeakCoupling::update(const unsigned long iterations)
    {
        for(unsigned long it=0l;it < iterations;it++){
            //TODO: if totalSign > threshold: break;
            steps 	+= 1;

            //auto sign	= 0;
            RealT t_n           = u(r_time);
            const RealT s_n     = static_cast<int>(u(r_spin) + 0.5);    	    	// auxiliary spin variable (-1 or 1)
            const RealT zeta 	= u(r_insert);
            const RealT zetap 	= u(r_accept);
            std::array<RealT,2> Sp;
            RealT A	= 0.0;

            //if((s_n != 0) and (s_n != 1)) LOG(WARNING) << "wrong aux spin value!";
            VLOG(2) << "Updating Configuration, current size: " << n << ", zetap: " << zetap;
            if (zeta < 0.5) {						    	                        // try to insert spin
                t_n *= config->beta;
                VLOG(3) << "Trying to add configuration, t_n : " << t_n << ", s_n: " << s_n;
                if(n == 0){						        	                // special treatment for empty M
                    //const RealT alpha = 0.5 + (2*s-1)*s_n*(0.5 + zeroShift);	                        // for 0 shift: \in {0,1}
                    // g0Call(0.0, s, s_n, 0.0)
                    
                    Sp[DOWN]  = g0Call(0.0, DOWN, s_n, 0.0);
                    Sp[UP]    = g0Call(0.0, UP, s_n, 0.0);
                    A = -Sp[DOWN]*Sp[UP]*config->beta*config->U;                                          // -b*U*detRatio/(n+1)
                    VLOG(3) << "Acceptance rate: " << A;
                    if(zetap < A){			                		                // propose insertion
                        VLOG(3) << "accepted";
                        M[UP] = Eigen::MatrixXd::Zero(1,1); M[DOWN] = Eigen::MatrixXd::Zero(1,1);
                        M[UP](0,0) = 1.0/Sp[UP];                M[DOWN](0,0) = 1.0/Sp[DOWN];
                        pushConfig(std::make_tuple(t_n, s_n));
                        currentSign = (2*(A>=0.)-1)*currentSign;
                        n = 1;
                        //sign = 2*(A>0)-1;
                    }
                } else {							    	                // n > 0, insertion
                    A = -config->beta*config->U/(static_cast<RealT>(n)+1);                                // -b*U*detRatio/(n+1)
                    RowVectorT R[2];
                    VectorT Q[2];
                    for(int s=0; s < _CONFIG_spins; s++){			                 	        // compute acceptance rate
                        RealT r_data[n];
                        RealT q_data[n];
                        for(int i=0; i<n; i++){
                            RealT tau = t_n - std::get<0>(confs[i]);
                            if(tau == 0.)
                            {
                                t_n += 0.00001;
                                tau += 0.00001;
                            }
                            //if(tau == 0.0) LOG(WARNING) << "off-diagonal G(0) sampled!";
                            r_data[i] = g0Call(t_n, s, s_n, std::get<0>(confs[i]));		        // (8.34) Generate new elements for M from Weiss Greens fct
                            q_data[i] = g0Call(std::get<0>(confs[i]), s, std::get<1>(confs[i]), t_n);
                        }
                        R[s] = Eigen::Map<RowVectorT>(r_data,n);
                        Q[s] = Eigen::Map<VectorT>(q_data,n);
                        RealT tmp     = g0Call(0., s, s_n, 0.) - R[s]*M[s]*Q[s];
                        Sp[s]= 1.0/tmp;	                                	// (8.39) 
                        A *= tmp;
                    }
                    VLOG(3) << "Acceptance rate: " << A;
                    if (zetap < A){ 		            	                	                // compute A(n+1 <- n) (8.36) (8.39)
                        VLOG(3) << "Accepted\n";
                        pushConfig(std::make_tuple(t_n, s_n));
                        for(int s=0; s < _CONFIG_spins; s++){			                 	        // compute acceptance rate
                            util::MInc(&(M[s]), R[s], Q[s], Sp[s], 0, false);
                            //sign = 2*(A>0)-1;
                        }
                        currentSign = (2*(A>=0.)-1)*currentSign;
                        n += 1;
                    }
                }										// inserted SConfig
                VLOG(4) << "After possible insertion";
            }
            else if ( n > 0 )
            {										// try to remove spin
                t_n *=n;
                VLOG(2) << "Trying to remove configuration";
                int rndConfPos = static_cast<int>(t_n);
                RealT A = -static_cast<RealT>(n)/(config->beta*config->U);
                Sp[DOWN] = M[DOWN](rndConfPos,rndConfPos);					// (8.37)
                Sp[UP] = M[UP](rndConfPos,rndConfPos);						// (8.37)
                A *= Sp[0]*Sp[1];
                if ( zetap < A){                    					//compute A(n-1 <- n) (8.37) (8.44) - WeakCoupling::acceptanceR
                    // TODO: loop over spins
                    VLOG(3) << "Accepted";
                    if(n == 1)
                    {
                        for(int s=DOWN; s < 2; s++){
                            M[s].resize(0, 0);
                        }
                        deleteConfig(0);
                        n = 0;
                    }
                    else
                    {
                        for(int s=DOWN; s < 2; s++){			                 	// compute acceptance rate
                            util::MDec(&(M[s]), rndConfPos);
                        }
                        deleteConfig(rndConfPos);
                        //sign = 2*(A>0)-1;
                        n -= 1;
                    }
                    currentSign = (2*(A>=0.)-1)*currentSign;
                }
            }
            if(A < 0)
                LOG(WARNING) << "Acceptance rate negative. A = " << A << ", n = " << n << ", k -> " << 2*(zeta>0.5)-1;
            if(n%50000 == 0)
            {
              for(int s=0;s<_CONFIG_spins;s++)
                  M[s] = rebuildM(s);
            }

            VLOG(3) << "before updateContribution. n = " << n;;
            if(_CONFIG_INT_DEBUG)
            {
              if(n>0){
              for(int s=0;s<_CONFIG_spins;s++){
                   auto Mr = rebuildM(s);
                   if(!Mr.isApprox(M[s], 1.e-8))
                   {
                       //LOG(ERROR) << "expected \n" << Mr << " \n but got \n " << M[s] << "\n difference: \n" << Mr - M[s] << "\n\n\n\n" <<  (Mr - M[s]).maxCoeff();
                       LOG(ERROR) <<"Fast matrix update failed in step " + std::to_string(steps);
                       //throw std::logic_error("Fast matrix update failed in step " + std::to_string(steps));
                   }
              }
            }
            }

            //updateContribution_OLD(1);
            updateContribution(1);


        }
    }


    /*! updates all variables for the MC simulation
     *  @param [in] sign +1/-1 or 0 to replicate the last sign (proposal not accepted => meassure last config again)
     */
    void WeakCoupling::updateContribution_OLD(int sign)
    {
        if(steps < burninSteps) return;	    // return while still in burn in period
        //TODO: tmp sign 
        if(!n)
        {
            std::vector<RealT> tmpVec = {static_cast<RealT>(sign)};
            config->world.send(0, static_cast<int>(MPI_MSG_TAGS::DATA), tmpVec);
            return;
        }
#ifdef MEASUREMENT_SHIFT
        const RealT rShift = config->beta*u(r_shift);
#else
        const RealT rShift = 0.0;
#endif
        VectorT out_d(n);
        VectorT out_u(n);
        std::vector<RealT> out(3*n+1);
        for(int i=0;i<n;i++)
        {
            out[i] = std::get<0>(confs[i]) - rShift;
            out_d(i) = g0Call(out[i],DOWN,std::get<1>(confs[i]),0.0);
            out_u(i) = g0Call(out[i],UP,std::get<1>(confs[i]),0.0);
        }
        //assume down == 0 
        Eigen::VectorXd::Map(&out[n], n) = (M[DOWN]*out_d).eval();
        Eigen::VectorXd::Map(&out[2*n], n) = (M[UP]*out_u).eval();
        out[3*n] = sign;
        //TODO: do local accumulator, later send sums over
        config->world.send(0, static_cast<int>(MPI_MSG_TAGS::DATA), out );
        return;
    }

    void WeakCoupling::updateContribution(int sign)
    {
        VLOG(3) << "entering update";
        if(steps < burninSteps) return;	    // return while still in burn in period

        currentSign = sign*currentSign;			    // remember last sign
        totalSign += 1;//currentSign;

        expOrdAcc(n, 0 );
        expOrdAcc(n, 1 );

        if(!n) return;
#ifdef MEASUREMENT_SHIFT
        const RealT zeta = config->beta*u(r_shift);
#else
        const RealT zeta = 0.0;
#endif

        VectorT tmpUP(n);
        VectorT tmpDOWN(n);
        for(int i=0;i<n;i++)
        {
            tmpDOWN(i) = g0Call(std::get<0>(confs[i]) - zeta,DOWN,std::get<1>(confs[i]), 0.0);
            tmpUP(i) = g0Call(std::get<0>(confs[i]) - zeta,UP,std::get<1>(confs[i]), 0.0);
        }
        tmpUP = M[UP]*tmpUP;
        tmpDOWN = M[DOWN]*tmpDOWN;
        for(int k=0;k<n;k++){		    // itBins(t)   = G0(t-t_k) * A_k
            RealT tau = std::get<0>(confs[k]) - zeta;
            const int sign2 = 2*(tau>=0.)-1;
            const int index = static_cast<int>(_CONFIG_maxSBins * (tau + (tau<0)*config->beta)/config->beta );
            itBinsDOWN[index](sign2*sign*tmpDOWN(k));
            itBinsUP[index](sign2*sign*tmpUP(k));
        }
        VLOG(3) << "leaving update";
    }

    MatrixT WeakCoupling::rebuildM(int spin)
    {
        //LOG(DEBUG) << "manually rebuilding M matrix";
        MatrixT Minv(n,n);
        for(int i=0;i<n;i++){
            //const RealT alpha = 0.5 + (2*spin-1)*std::get<1>(confs[i])*(0.5 + zeroShift);
            Minv(i,i) = g0Call(0.0, spin, std::get<1>(confs[i]),0.0);
            for(int j=i+1;j<n;j++){
                Minv(i,j) =  g0Call(std::get<0>(confs[i]), spin, std::get<1>(confs[i]), std::get<0>(confs[j])) ;
                Minv(j,i) =  g0Call(std::get<0>(confs[j]), spin, std::get<1>(confs[j]), std::get<0>(confs[i])) ;
            }
        }
        return Minv.inverse();
    }

} //end namspace DMFT
