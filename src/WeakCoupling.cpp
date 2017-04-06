#include "WeakCoupling.hpp"
//TODO: use BOOST_FOR and BOOST_MPI
namespace DMFT
{

    WeakCoupling::WeakCoupling(
            GreensFct &g0, GreensFct &gImp, const Config& config, const RealT zeroShift, const unsigned int burninSteps
            ):
        config(config), g0(g0), gImp(gImp), zeroShift(zeroShift), burninSteps(burninSteps) 
    {
        n 		= 0;
        steps 		= 0;
        totalSign	= 0;
        lastSign	= 1;
        itBins		= ImTG::Zero(_CONFIG_maxSBins*config.beta, 2);
        r_time.split(5, 0);
        r_spin.split(5, 1);
        r_insert.split(5, 2);
        r_accept.split(5, 3);
        r_shift.split(5, 4);
        r_time.split(config.local.size() , config.local.rank());           // choose sub−stream no. rank out of size streams
        r_spin.split(config.local.size() , config.local.rank());           // choose sub−stream no. rank out of size streams
        r_insert.split(config.local.size() , config.local.rank());         // choose sub−stream no. rank out of size streams
        r_accept.split(config.local.size() , config.local.rank());         // choose sub−stream no. rank out of size streams
        r_shift.split(config.local.size() , config.local.rank());         // choose sub−stream no. rank out of size streams
    }


    WeakCoupling::~WeakCoupling()
    {
    }


    void WeakCoupling::computeImTGF(void){
#ifdef WK_NAIVE_MEASUREMENT
        gImp.setByT(itBins2);
#else
        for(int s=0;s<2;s++){
            for(RealT t=0.0;t<config.beta;t+=1.0/_CONFIG_maxTBins){
                RealT sumIt = 0.0;
                //TODO: to stable boost accumulation here
                for(int j=0;j<_CONFIG_maxSBins*config.beta;j++){
                    if(itBins(j) == 0.0) continue;
                    const RealT tp = static_cast<RealT>(j)/_CONFIG_maxSBins;
                    const RealT gfVal = g0(t-tp,s);
                    sumIt += gfVal*itBins(j);
                }
                gImp.setByT(t,s, g0(t,s) - sumIt/totalSign);
            }
        }
#endif
    }

    void WeakCoupling::computeMatGreensFct(void)
    {
        //TODO: flag for meassurement of mats G
        for(unsigned int n=0;n<_CONFIG_maxMatsFreq; n++){
            ComplexT sumWn = 0.0;
            RealT mfreq = mFreq(n,config.beta);
            for(int s=0;s<2;s++){
                for(int j=0;j<_CONFIG_maxSBins*config.beta;j++){
                    if(itBins(j,s) == 0.0) continue;
                    const RealT tp = static_cast<RealT>(j)/_CONFIG_maxSBins;
                    sumWn += std::exp(ComplexT(0.0, mfreq*tp))*itBins(j,s);
                }
                gImp.setByMFreq(n,s, g0.getByMFreq(n,s) - g0.getByMFreq(n,s)*sumWn/static_cast<RealT>(totalSign));
            }
        }
    }

    //TODO:	do this properly: accumulator, overflow, kahan, SBin not *beta, mats != tbins 
    void WeakCoupling::computeImpGF(void)
    {
        (config.local.barrier)();                                           // wait for all MC runners to finish
        if(config.local.rank() == 0)
        config.world.send(0, static_cast<int>(MPI_MSG_TAGS::COMM_END));     // inform accumulators, that we are done

        for(int s=0;s<2;s++){
            for(int n=0;n<_CONFIG_maxMatsFreq; n++){
                RealT t 		= config.beta*static_cast<RealT>(n)/_CONFIG_maxMatsFreq;
                ComplexT sumWn	= 0.0;
                RealT sumIt		= 0.0;
                RealT mfreq		= mFreq(n,config.beta);
                for(int j=0;j<_CONFIG_maxSBins*config.beta;j++){
                    const RealT bVal = itBins(j,s);
                    if(bVal == 0.0) continue;
                    const RealT tp = static_cast<RealT>(j)/_CONFIG_maxSBins;
                    sumWn += std::exp(ComplexT(0.0, mfreq*tp))*bVal;
                    sumIt += g0(t-tp,s)*bVal;
                }
                gImp.setByMFreq(n,s, g0.getByMFreq(n,s) - g0.getByMFreq(n,s)*sumWn/static_cast<RealT>(totalSign));
                gImp.setByT(t,s, g0(t,s) - sumIt/totalSign);
            }
        }
    }


    void WeakCoupling::update(const unsigned long iterations)
    {
        for(unsigned long it=0l;it < iterations;it++){
            steps 	+= 1;

            //auto sign	= 0;
            RealT t_n     = u(r_time);
            const RealT s_n     = 2*static_cast<int>(u(r_spin) + 0.5) - 1;    	    	// auxiliary spin variable (-1 or 1)
            const RealT zeta 	= u(r_insert);
            const RealT zetap 	= u(r_accept);
            std::array<RealT,2> Sp;
            RealT A	= 0.0;

            VLOG(2) << "Updating Configuration, current size: " << n << ", zetap: " << zetap;
            if (zeta < 0.5) {						    	// try to insert spin
                VLOG(3) << "Trying to add configuration";
                t_n *= config.beta;
                if(n == 0){						        	// special treatment for empty M
                    //const RealT alpha = 0.5 + (2*s-1)*s_n*(0.5 + zeroShift);	// for 0 shift: \in {0,1}
                    //TODO: workaround (performance): assume SPIN::DOWN == 0
                    Sp[DOWN]  = g0(0.0,DOWN) - 0.5 + s_n*(0.5 + zeroShift);
                    Sp[UP]    = g0(0.0,UP) - 0.5 - s_n*(0.5 + zeroShift);
                    A = -Sp[0]*Sp[1]*config.beta*config.U;                              // -b*U*detRatio/(n+1)
                    VLOG(3) << "Acceptance rate: " << A;
                    if(zetap < A){			                		// propose insertion
                        VLOG(3) << "accepted";
                        M[UP] = Eigen::MatrixXd::Zero(1,1); M[DOWN] = Eigen::MatrixXd::Zero(1,1);
                        M[UP] << 1.0/Sp[UP];                M[DOWN] << 1.0/Sp[DOWN];
                        pushConfig(std::make_tuple(t_n, s_n));
                        n = 1;
                        //sign = 2*(A>0)-1;
                    }
                } else {							    	// n > 0, insertion
                    A = -config.beta*config.U/(static_cast<RealT>(n)+1);                // -b*U*detRatio/(n+1)
                    RowVectorT R[2];
                    VectorT Q[2];
                    RealT tmp;

                    // assume down == 0 for performance
                    for(int s=DOWN; s < 2; s++){			                 	// compute acceptance rate
                        RealT r_data[n];
                        RealT q_data[n];
                        for(int i=0; i<n;i++ ){
                            r_data[i] = g0(t_n- std::get<0>(confs[i]),s);		// (8.34) Generate new elements for M from Weiss Greens fct
                            q_data[i] = g0(std::get<0>(confs[i])-t_n,s);
                        }
                        R[s] = Eigen::Map<RowVectorT>(r_data,n);
                        Q[s] = Eigen::Map<VectorT>(q_data,n);
                        //const RealT alpha = 0.5 + (2*s-1)*s_n*(0.5 + zeroShift);	// for 0 shift: \in {0,1}
                        RealT tmp     = (g0(0.0,s) - 0.5 + (2*s-1)*s_n*(0.5 + zeroShift)) - R[s]*M[s]*Q[s];
                        Sp[s]= 1.0/tmp;	                                	// (8.39) 
                        A *= tmp;
                    }

                    VLOG(3) << "Acceptance rate: " << A;
                    if (zetap < A){ 		            	                	// compute A(n+1 <- n) (8.36) (8.39)
                        VLOG(3) << "Accepted\n";
                        pushConfig(std::make_tuple(t_n, s_n));
                        // assume down == 0 for performance
                        for(int s=DOWN; s < 2; s++){			                 	// compute acceptance rate
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
                RealT A = -static_cast<RealT>(n)/(config.beta*config.U);   //QUESTION: - missing here, go through proof again
                Sp[DOWN] = M[DOWN](rndConfPos,rndConfPos);					// (8.37)
                Sp[UP] = M[UP](rndConfPos,rndConfPos);						// (8.37)
                A *= Sp[0]*Sp[1];
                if ( zetap < A){                    					//compute A(n-1 <- n) (8.37) (8.44) - WeakCoupling::acceptanceR
                    // TODO: loop over spins
                    VLOG(3) << "Accepted";

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

            VLOG(3) << "before updateContribution. n = " << n;
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

            updateContribution(1); //TODO: for now assume no sign problem
        }
    }


    /*! updates all variables for the MC simulation
     *  @param [in] sign +1/-1 or 0 to replicate the last sign (proposal not accepted => meassure last config again)
     */
    void WeakCoupling::updateContribution(int sign)
    {
        if(steps < burninSteps) return;	    // return while still in burn in period
        if(!n)
        {
            config.world.send(0, static_cast<int>(MPI_MSG_TAGS::DATA), sign );
            return;
        }
#ifdef MEASUREMENT_SHIFT
        const RealT rShift = config.beta*u(r_shift);
#else
        const RealT rShift = 0.0;
#endif
        MatrixT tmp(n,3);
        //TODO: OPTIMIZE: this is horribly slow
        for(int i=0;i<n;i++)
        {
            tmp(i,0) = std::get<0>(confs[i]) - rShift;
            tmp(i,1) = g0(tmp(i,0),UP);
            tmp(i,2) = g0(tmp(i,0),DOWN);
        }
        tmp.col(1) = (M[UP]*tmp.col(1)).eval();
        tmp.col(2) = (M[DOWN]*tmp.col(2)).eval();
        std::vector<RealT> tmpVec(tmp.data(), tmp.data() + tmp.rows()*tmp.cols());
        tmpVec.push_back(static_cast<RealT>(sign));

        //TODO: overload boost serialize to directly send eigen arrays
        config.world.send(0, static_cast<int>(MPI_MSG_TAGS::DATA), tmpVec );
        return;
    }

    MatrixT WeakCoupling::rebuildM(int spin)
    {
        LOG(DEBUG) << "manually rebuilding M matrix";
        MatrixT Minv(n,n);
        for(int i=0;i<n;i++){
            const RealT alpha = 0.5 + (2*spin-1)*std::get<1>(confs[i])*(0.5 + zeroShift);
            Minv(i,i) = g0(0.0, spin) - alpha;
            for(int j=i+1;j<n;j++){
                Minv(i,j) =  g0(std::get<0>(confs[i]) - std::get<0>(confs[j]), spin) ;
                Minv(j,i) =  g0(std::get<0>(confs[j]) - std::get<0>(confs[i]), spin) ;
            }
        }
        return Minv.inverse();
    }

}
