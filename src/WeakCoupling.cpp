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
        lastSign	= 1;
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
        for(long j = 0; j< binCount; j++)
        {
            RealT tp = j*config->beta/(binCount);
            const RealT bValUP = boost::accumulators::sum(itBinsUP[j]);
            const RealT bValDOWN = boost::accumulators::sum(itBinsDOWN[j]);
            if(bValUP)
            {
                for(long i=0; i<itCount; i++)
                {
                    RealT t = i*config->beta/(itCount);
                    g_it(i,UP) += (*g0)(t-tp,UP)*bValUP;
                }
#ifdef MATSUBARA_MEASUREMENT
                for(long k=0; k<mfCount;k++)
                {
                    RealT mfreq = gImp->isSymmetric() ? mFreqS(k,config->beta) : mFreq(k, config->beta);
                    g_wn(k,UP) += std::exp(ComplexT(0.0,mfreq*tp))*bValUP;
                }
#endif
            }
            if(bValDOWN)
            {
                for(long i=0; i<itCount; i++)
                {
                    RealT t = i*config->beta/(itCount);
                    g_it(i,DOWN) += (*g0)(t-tp,DOWN)*bValDOWN;
                }
#ifdef MATSUBARA_MEASUREMENT
                for(long k=0; k<mfCount;k++)
                {
                    RealT mfreq = gImp->isSymmetric() ? mFreqS(k,config->beta) : mFreq(k, config->beta);
                    g_wn(k,DOWN) += std::exp(ComplexT(0.0,mfreq*tp))*bValDOWN;
                }
#endif
            }
        }
        g_it = g0->getItGF() - g_it/totalSign;
        gImp->setByT(g_it);
#ifndef MATSUBARA_MEASUREMENT
        gImp->transformTtoM();
#else
        g_wn = g0->getMGF() - g0->getMGF()*g_wn/totalSign;
        gImp->setByMFreq(g_wn);
        //gImp->transformMtoT();
#endif
        //gImp->transformMtoT();
        reset();
    }

    void WeakCoupling::reset()
    {
        n = 0;
        totalSign = 0;
        lastSign = 1;
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
        gfCache[UP] = tmp2;
        gfCache[DOWN] = tmp2;
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
                    
                    Sp[DOWN]  = g0Call(0.0, s_n, DOWN, 0.0);
                    Sp[UP]    = g0Call(0.0, s_n, UP, 0.0);
                    A = -Sp[DOWN]*Sp[UP]*config->beta*config->U;                                          // -b*U*detRatio/(n+1)
                    if(A < 0 && A < -0.01) LOG(WARNING) << "Acceptance rate for insertion negative. n = 0, A=" << A; //Sp[DOWN]: " << Sp[DOWN] << ", SP[UP]: " << Sp[UP];
                    VLOG(3) << "Acceptance rate: " << A;
                    if(zetap < A){			                		                // propose insertion
                        VLOG(3) << "accepted";
                        M[UP] = Eigen::MatrixXd::Zero(1,1); M[DOWN] = Eigen::MatrixXd::Zero(1,1);
                        M[UP] << 1.0/Sp[UP];                M[DOWN] << 1.0/Sp[DOWN];
                        pushConfig(std::make_tuple(t_n, s_n));
                        n = 1;
                        //sign = 2*(A>0)-1;
                    }
                } else {							    	                // n > 0, insertion
                    A = -config->beta*config->U/(static_cast<RealT>(n)+1);                                // -b*U*detRatio/(n+1)
                    RowVectorT R[2];
                    VectorT Q[2];
                    RealT tmp;

                    // assume down == 0 for performance
                    for(int s=DOWN; s < 2; s++){			                 	        // compute acceptance rate
                        RealT r_data[n];
                        RealT q_data[n];
                        for(int i=0; i<n;i++ ){
                            const RealT tau = t_n - std::get<0>(confs[i]);
                            //if(tau == 0.0) LOG(WARNING) << "off-diagonal G(0) sampled!";
                            r_data[i] = g0Call_od(t_n, s, std::get<0>(confs[i]));		        // (8.34) Generate new elements for M from Weiss Greens fct
                            q_data[i] = g0Call_od(std::get<0>(confs[i]), s,t_n);
                        }
                        R[s] = Eigen::Map<RowVectorT>(r_data,n);
                        Q[s] = Eigen::Map<VectorT>(q_data,n);
                        RealT tmp     = g0Call(0., s, s_n, 0.) - R[s]*M[s]*Q[s];
                        Sp[s]= 1.0/tmp;	                                	// (8.39) 
                        A *= tmp;
                    }
                    //if(A<0 && A < -0.01) LOG(WARNING) << "Acceptance rate for insertion negative. n = " << n << ", A = " << A;
                    /*  if(A<0) LOG(WARNING) << "Acceptance rate for insertion negative. n = " << n << ", A = " << A << ", s_n: " << s_n
                        << "\n (*g0)(0, DOWN): " << (*g0)(0.0,DOWN) << " -> " << (*g0)(0.0,DOWN) - (0.5 + (2*(s_n==DOWN)-1)*zeroShift) <<  ", tmp[down]: " << R[DOWN]*M[DOWN]*Q[DOWN]
                        << "\n  (*g0)(0, UP): " << (*g0)(0.0,UP)  << " -> " << (*g0)(0.0,UP) - (0.5 + (2*(s_n==UP)-1)*zeroShift) << ", tmp[up]: " << R[UP]*M[UP]*Q[UP];
                        */
                    VLOG(3) << "Acceptance rate: " << A;
                    if (zetap < A){ 		            	                	                // compute A(n+1 <- n) (8.36) (8.39)
                        VLOG(3) << "Accepted\n";
                        pushConfig(std::make_tuple(t_n, s_n));
                        // assume down == 0 for performance
                        for(int s=DOWN; s < 2; s++){			                 	        // compute acceptance rate
                            M[s].conservativeResize(n+1,n+1);
                            M[s](n,n)			         =  Sp[s];               				// (8.40)
                            M[s].topRightCorner(n,1).noalias()	 = -(M[s].topLeftCorner(n,n)*Q[s])*Sp[s];		// (8.41)
                            M[s].bottomLeftCorner(1,n).noalias() = -Sp[s]*(R[s]*M[s].topLeftCorner(n,n));
                            M[s].topLeftCorner(n,n)		 += M[s].topRightCorner(n,1) * M[s].bottomLeftCorner(1,n)/Sp[s];
                            //sign = 2*(A>0)-1;
                        }
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
                RealT A = -static_cast<RealT>(n)/(config->beta*config->U);   //QUESTION: - missing here, go through proof again
                Sp[DOWN] = M[DOWN](rndConfPos,rndConfPos);					// (8.37)
                Sp[UP] = M[UP](rndConfPos,rndConfPos);						// (8.37)
                A *= Sp[0]*Sp[1];
                if(A<0 && A < -0.01) LOG(WARNING) << "Acceptance rate for insertion negative. n = " << n << ", A = " << A;
                if ( zetap < A){                    					//compute A(n-1 <- n) (8.37) (8.44) - WeakCoupling::acceptanceR
                    // TODO: loop over spins
                    VLOG(3) << "Accepted";
                    if(n == 1)
                    {
                        MatrixT tmp;
                        VectorT tmp2;
                        for(int s=DOWN; s < 2; s++){
                            M[s] = tmp;
                            gfCache[s] = tmp2;
                        }
                        n = 0;
                    }
                    else
                    {
                        for(int s=DOWN; s < 2; s++){			                 	// compute acceptance rate
                            if(rndConfPos < n-1){
                                M[s].row(rndConfPos).swap(M[s].row(n-1));
                                M[s].col(rndConfPos).swap(M[s].col(n-1));
                            }
                            M[s].topLeftCorner(n-1,n-1) = (M[s].topLeftCorner(n-1,n-1) - M[s].topRightCorner(n-1,1)
                                    * M[s].bottomLeftCorner(1,n-1) / M[s](n-1,n-1)).eval();	//M = P-Q*R/S (8.45)
                            M[s].conservativeResize(n-1,n-1);
                        }
                        swapConfigs(rndConfPos, n-1);    								// swap the same configs as in Ms
                        popConfig();			   									// pop last one 
                        //sign = 2*(A>0)-1;
                        n -= 1;
                    }
                }
            }

            VLOG(3) << "before updateContribution. n = " << n;;
            /*#ifdef DEBUG_MODE
              if(n>0){
              acceptanceR();
              for(int s=0;s<_CONFIG_spins;s++){
              LOG(DEBUG) << "rebuilt M: " << rebuildM(s);
              LOG(DEBUG) << "actual M: " << Ms[s];
              }
            //static_assert( (Ms[0] - rebuildM(0)).norm() < 0.001, "Ms[0] matrix not equal to explicit inverse");
            //static_assert( (Ms[1] - rebuildM(1)).norm() < 0.001, "Ms[1] matrix not equal to explicit inverse");
            }
#endif*/

            //updateContribution_OLD(1); //TODO: for now assume no sign problem
            updateContribution(1); //TODO: for now assume no sign problem


        }
    }


    /*! updates all variables for the MC simulation
     *  @param [in] sign +1/-1 or 0 to replicate the last sign (proposal not accepted => meassure last config again)
     */
    void WeakCoupling::updateContribution_OLD(int sign)
    {
        if(steps < burninSteps) return;	    // return while still in burn in period
        //TODO: tmp sign 
        if(!sign) sign = lastSign;
        lastSign = sign;
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

        //if(!sign) sign = lastSign;		    // update got rejected, use last sign
        //lastSign = sign;			    // remember last sign
        totalSign += 1;//sign;
        expOrdAcc(M[0].rows(), 0 );
        expOrdAcc(M[1].rows(), 1 );

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
            const int sign2 = 2*(tau>0)-1;
            const int index = static_cast<int>(_CONFIG_maxSBins * (tau + (tau<0)*config->beta)/config->beta );
            itBinsDOWN[index](sign2*sign*tmpDOWN(k));
            itBinsUP[index](sign2*sign*tmpUP(k));
        }
        VLOG(3) << "leaving update";
    }

    MatrixT WeakCoupling::rebuildM(int spin)
    {
        LOG(DEBUG) << "manually rebuilding M matrix";
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
