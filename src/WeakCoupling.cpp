#include "WeakCoupling.hpp"
//TODO: use BOOST_FOR and BOOST_MPI
namespace DMFT
{

    WeakCoupling::WeakCoupling(GreensFct &g0, GreensFct &gImp, const Config& config, const RealT zeroShift, const unsigned int burninSteps):
        config(config), g0(g0), gImp(gImp), zeroShift(zeroShift), burninSteps(burninSteps), MPI_size(MPI::COMM_WORLD.Get_size()), MPI_rank(MPI::COMM_WORLD.Get_rank())
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
        r_time.split(MPI_size, MPI_rank);           // choose sub−stream no. rank out of size streams
        r_spin.split(MPI_size, MPI_rank);           // choose sub−stream no. rank out of size streams
        r_insert.split(MPI_size, MPI_rank);         // choose sub−stream no. rank out of size streams
        r_accept.split(MPI_size, MPI_rank);         // choose sub−stream no. rank out of size streams
        r_shift.split(MPI_size, MPI_rank);         // choose sub−stream no. rank out of size streams
    }


    WeakCoupling::~WeakCoupling()
    {
    }

    void WeakCoupling::acceptanceR()
    {
        RealT deleteRate = -(n+1)/(config.beta*config.U);
        RealT insertRate = -config.beta*config.U/(n+1);
        RowVectorT R[2];
        VectorT Q[2];
        for(int pos=0;pos<n;pos++){
            for(int s=0; s < 2; s++){
                RealT r_data[n];
                RealT q_data[n];

                for(int i=0; i<n;i++ ){
                    r_data[i] = g0(std::get<0>(confs[pos]) - std::get<0>(confs[i]),s);
                    q_data[i] = g0(std::get<0>(confs[i]) - std::get<0>(confs[pos]),s);
                }
                R[s] = Eigen::Map<RowVectorT>(r_data,n);
                Q[s] = Eigen::Map<VectorT>(q_data,n);
            }
            const RealT alphaDown   = 0.5 + std::get<1>(confs[pos])*(0.5 + zeroShift);
            const RealT alphaUp     = 0.5 - std::get<1>(confs[pos])*(0.5 + zeroShift);
            insertRate *= ((g0(0.0,0) - alphaDown) - R[0]*Mdown*Q[0]);
            insertRate *= ((g0(0.0,1) - alphaUp) - R[1]*Mup*Q[1]);

            deleteRate *= Mup(pos,pos);
            deleteRate *= Mdown(pos,pos);
            LOG(DEBUG) << "n: "<<n<<", insertion acceptance: " << insertRate << ", removal acceptance: " << deleteRate << ", ratio: "<< insertRate/deleteRate;
        }
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
        for(int s=0;s<_CONFIG_spins;s++){
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
                //TODO: better accumulator
                //if(std::isnan(std::real(sumWn)) or std::isinf(std::real(sumWn)) or std::isnan(std::imag(sumWn)) or std::isinf(std::imag(sumWn)) ) LOG(ERROR) << "Overflow during computation of G_Imp(i wn)";
                //if(std::isnan(sumIt) or std::isinf(sumIt)) LOG(ERROR) << "Overflow during computation of G_Imp(tau)";
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
            RealT 	SpUp,SpDown;
            RealT A	= 0.0;

            VLOG(2) << "Updating Configuration, current size: " << n << ", zetap: " << zetap;
            if (zeta < 0.5) {						    	// try to insert spin
                VLOG(3) << "Trying to add configuration";
                t_n *= config.beta;
                if(n == 0){						        	// special treatment for empty M
                    //const RealT alpha = 0.5 + (2*s-1)*s_n*(0.5 + zeroShift);	// for 0 shift: \in {0,1}
                    SpDown  = g0(0.0,0) - 0.5 - s_n*(0.5 + zeroShift);
                    SpUp    = g0(0.0,1) - 0.5 + s_n*(0.5 + zeroShift);                  // s==1 => spin down
                    A = -SpUp*SpDown*config.beta*config.U;                              // -b*U*detRatio/(n+1)
                    VLOG(3) << "Acceptance rate: " << A;
                    if(zetap < A){			                		// propose insertion
                        VLOG(3) << "accepted";
                        Mup = Eigen::MatrixXd::Zero(1,1);
                        Mup << 1.0/SpUp;
                        Mdown = Eigen::MatrixXd::Zero(1,1);
                        Mdown << 1.0/SpDown;
                        pushConfig(std::make_tuple(t_n, s_n));
                        n = 1;
                        //sign = 2*(A>0)-1;
                    }
                } else {							    	// n > 0, insertion
                    A = -config.beta*config.U/(static_cast<RealT>(n)+1);                // -b*U*detRatio/(n+1)
                    RowVectorT R[2];
                    VectorT Q[2];

                    for(int s=0; s < 2; s++){			                 	// compute acceptance rate
                        RealT r_data[n];
                        RealT q_data[n];
                        //const RealT alpha = 0.5 + (2*s-1)*s_n*(0.5 + zeroShift);	// for 0 shift: \in {0,1}

                        for(int i=0; i<n;i++ ){
                            r_data[i] = g0(t_n- std::get<0>(confs[i]),s);		// (8.34) Generate new elements for M from Weiss Greens fct
                            q_data[i] = g0(std::get<0>(confs[i])-t_n,s);
                        }
                        R[s] = Eigen::Map<RowVectorT>(r_data,n);
                        Q[s] = Eigen::Map<VectorT>(q_data,n);
                    }
                    RealT tmpDown   = (g0(0.0,0) - 0.5 - s_n*(0.5 + zeroShift)) - R[0]*Mdown*Q[0];
                    RealT tmpUp     = (g0(0.0,1) - 0.5 + s_n*(0.5 + zeroShift)) - R[1]*Mup*Q[1];
                    SpDown= 1.0/tmpDown;	                                	// (8.39) 
                    SpUp= 1.0/tmpUp;	                                         	// (8.39) 
                    A *= tmpUp*tmpDown;


                    VLOG(3) << "Acceptance rate: " << A;
                    if (zetap < A){ 		            	                	// compute A(n+1 <- n) (8.36) (8.39)
                        VLOG(3) << "Accepted\n";
                        pushConfig(std::make_tuple(t_n, s_n));
                        Mup.conservativeResize(n+1,n+1);
                        Mdown.conservativeResize(n+1,n+1);
                        Mdown(n,n)			        = SpDown;               				// (8.40)
                        Mdown.topRightCorner(n,1).noalias()	= -(Mdown.topLeftCorner(n,n)*Q[0])*SpDown;		// (8.41)
                        Mdown.bottomLeftCorner(1,n).noalias()   = -SpDown*(R[0]*Mdown.topLeftCorner(n,n));
                        Mdown.topLeftCorner(n,n)		+= Mdown.topRightCorner(n,1) * Mdown.bottomLeftCorner(1,n)/SpDown;
                        Mup(n,n)			        = SpUp;			                        	// (8.40)
                        Mup.topRightCorner(n,1).noalias()	= -(Mup.topLeftCorner(n,n)*Q[1])*SpUp;  		// (8.41)
                        Mup.bottomLeftCorner(1,n).noalias()     = -SpUp*(R[1]*Mup.topLeftCorner(n,n));
                        Mup.topLeftCorner(n,n)		        += Mup.topRightCorner(n,1) * Mup.bottomLeftCorner(1,n)/SpUp;
                        //sign = 2*(A>0)-1;
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
                    SpDown = Mdown(rndConfPos,rndConfPos);					// (8.37)
                    SpUp = Mup(rndConfPos,rndConfPos);						// (8.37)
                    A *= SpUp*SpDown;
                    if ( zetap < A){                    					//compute A(n-1 <- n) (8.37) (8.44) - WeakCoupling::acceptanceR
                        VLOG(3) << "Accepted";
                        if(rndConfPos < n-1){
                            Mdown.row(rndConfPos).swap(Mdown.row(n-1));
                            Mdown.col(rndConfPos).swap(Mdown.col(n-1));

                            Mup.row(rndConfPos).swap(Mup.row(n-1));
                            Mup.col(rndConfPos).swap(Mup.col(n-1));
                        }
                        Mup.topLeftCorner(n-1,n-1) = (Mup.topLeftCorner(n-1,n-1) - Mup.topRightCorner(n-1,1)
                                * Mup.bottomLeftCorner(1,n-1) / Mup(n-1,n-1)).eval();	//M = P-Q*R/S (8.45)
                        Mup.conservativeResize(n-1,n-1);
                        Mdown = (Mdown.topLeftCorner(n-1,n-1) - Mdown.topRightCorner(n-1,1)
                                * Mdown.bottomLeftCorner(1,n-1) / Mdown(n-1,n-1)).eval();	//M = P-Q*R/S (8.45)
                        Mdown.conservativeResize(n-1,n-1);
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
            VLOG(3) << "entering update";
            if(steps < burninSteps) return;	    // return while still in burn in period

            //if(!sign) sign = lastSign;		    // update got rejected, use last sign
            //lastSign = sign;			    // remember last sign
            totalSign += 1;//sign;


            if(!n) return;
#ifdef MEASUREMENT_SHIFT
            const RealT zeta = config.beta*u(r_shift);
#else
            const RealT zeta = 0.0;
#endif

            VectorT gTmpUP(n);
            VectorT gTmpDOWN(n);
            for(int i=0;i<n;i++)
            {
                gTmpUP(i) = g0(std::get<0>(confs[i]) - zeta,0);
                gTmpDOWN(i) = g0(std::get<0>(confs[i]) - zeta,1);
            }
#ifdef MEASUREMENT_SHIFT 
            VectorT tmpUP = Mup*gTmpUP;
            VectorT tmpDOWN = Mdown*gTmpDOWN;
#else
            VectorT tmpUP = Mup*gTmpUP;
            VectorT tmpDOWN = Mdown*gTmpDOWN;
#endif
#ifdef WK_NAIVE_MEASUREMENT
            for(auto t=0.0;t<config.beta;t+=config.beta/_CONFIG_maxTBins)
            {
                for(int k=0;k<n;k++)
                {
                    const auto tau = std::get<0>(confs[k]);
                    itBins2(t,0) += g0(t-tau, 0)*tmpUP(k);
                    itBins2(t,1) += g0(t-tau, 1)*tmpDOWN(k);
                }
            }
#else
            //TODO: accumulator with arbitrary precision
            for(int k=0;k<n;k++){		    // itBins(t)   = G0(t-t_k) * A_k
                RealT tau = std::get<0>(confs[k]) - zeta;
                const int sign2 = 2*(tau>0)-1;
                itBins(static_cast<int>(_CONFIG_maxSBins * (tau + (tau<0)*config.beta) ), 0) += sign2*sign*tmpDOWN(k);
                itBins(static_cast<int>(_CONFIG_maxSBins * (tau + (tau<0)*config.beta) ), 1) += sign2*sign*tmpUP(k);
            }
#endif
            VLOG(3) << "leaving update";
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
