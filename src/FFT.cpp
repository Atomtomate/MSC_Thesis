#include "FFT.hpp"
#ifndef IOHELPER_H_
#include "IOhelper.hpp"
#endif
namespace DMFT
{

        //TODO: this can all be obtained from GrennFct object directly
        ComplexT FFT::tail(const RealT wn, const std::vector<std::array<RealT,2> >& tail, const int spin) const
        {
            ComplexT iwn(0.0, wn);
            //for(int i = 0; i < tail.size(); i++)
            //{
            //    res += tail[i][spin]*(std::pow(iwn, -i));
            //}
            ComplexT res = tail[0][spin] + tail[1][spin]/iwn - tail[2][spin]/(wn*wn) - tail[3][spin]/(iwn*wn*wn) + tail[4][spin]/(wn*wn*wn*wn);
            //return 1./iwn;
            return res;
        }

        //TODO: this can all be obtained from GrennFct object directly
        RealT FFT::ftTail(const RealT tau, const std::vector<std::array<RealT,2> >& tail, const int spin) const
        {
            //return -0.5;
            VLOG(8) << "ft Tail correction: " <<-tail[1][spin]/2.0 + tail[2][spin]*(2.0*tau - _beta)/4.0 + tail[3][spin]*tau*(_beta-tau)/4.0\
                + tail[4][spin]*(2.0*tau - _beta)*(2.0*tau*tau - 2.0*tau*_beta - _beta*_beta)/48.0;
            return tail[0][spin] - tail[1][spin]/2.0 + tail[2][spin]*(2.0*tau - _beta)/4.0 + tail[3][spin]*tau*(_beta-tau)/4.0\
                + tail[4][spin]*(2.0*tau - _beta)*(2.0*tau*tau - 2.0*tau*_beta - _beta*_beta)/48.0;
        }


        void FFT::transformMtoT(const MatG &from, ImTG &to, bool symmetric)
        {
            LOG(WARNING) << "untested";
            std::array<RealT,2> z = {.0, .0};
            //TODO: # tail coeffs in Config?
            const std::vector<std::array<RealT,2> > t = {z,z,z,z};
            transformMtoT(from, to, t, symmetric);
        }

        //TODO: use http://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html
        void FFT::transformMtoT(const MatG &from, ImTG &to, const std::vector<std::array<RealT,2> >& tail_coeff, bool symmetric)
        {
            //TODO: this is a stupid workaround because the plan creation is not thread safe
            const int tBins = _CONFIG_maxTBins;
            fftw_complex fftw3_input[tBins];
            fftw_complex fftw3_output[tBins];
            auto planMtoT = fftw_plan_dft_1d(tBins,fftw3_input,fftw3_output,FFTW_FORWARD,FFTW_ESTIMATE);
            VLOG(5) << "transforming M to T, coeffs: " << tail_coeff[0][0] << ", " << tail_coeff[1][0] << ", " << tail_coeff[2][0] << ", " << tail_coeff[3][0] << ", " << tail_coeff[4][0] << ", " << tail_coeff[5][0];
            // setting parameters
            const RealT loc_fft_dt   = _beta/(tBins); // MARK: + 1
            const RealT loc_fft_tmin =  0.;//0.5*_beta/tBins;
            const int nMF = from.rows();
            const int n_min = 0;
            const int n_max = n_min + nMF;
            const RealT n_fac = symmetric ? 2 : 1;
            const RealT loc_fft_wmin = mFreqS( n_min, _beta);
            for(int s=0;s<2;s++){
                if(std::abs(tail_coeff[1][s] - 1.0) > 0.1) LOG(WARNING) << "unexpected tail coefficient, only symmetric GF tested! Expected 1, got " << tail_coeff[1][s];
                // preparing fftw input
                memset(fftw3_input, 0, tBins*sizeof(fftw_complex));
                for(int n = n_min; n < n_max; n++)          // construct negative wn for symmetric GF
                {
                    const RealT wn = mFreqS(n, _beta);
                    ComplexT input_tmp =  (n_fac*std::exp( ComplexT(0.0,-loc_fft_tmin * wn ))*(from(n+n_min,s) + 1.0*ComplexT(0.,tail_coeff[1][s]/wn)))/_beta; // - tail_coeff[1][s]/ComplexT(0.,wn) 
                    fftw3_input[n+n_min][0] = input_tmp.real();
                    fftw3_input[n+n_min][1] = input_tmp.imag();
                }

                // transforming
                fftw_execute(planMtoT);

                // preparing output
                if(false)
                {
                    for(int t=0; t<tBins/2; t++){
                        const RealT tau_k = loc_fft_dt*(t); //MARK: + 1
                        ComplexT tmp = ComplexT(fftw3_output[t][0],fftw3_output[t][1]);
                        to(t,s) =  (tmp*std::exp(ComplexT(0.0, -loc_fft_wmin*(tau_k - loc_fft_tmin)))).real() - 0.5*tail_coeff[1][s]; 
                        to(to.rows()-t-1,s) = to(t,s);
                    }
                }else{
                    for(int t=0; t<tBins; t++){
                        const RealT tau_k = loc_fft_dt*(t); //MARK: + 1
                        ComplexT tmp = ComplexT(fftw3_output[t][0],fftw3_output[t][1]);
                        to(t,s) =  (tmp*std::exp(ComplexT(0.0, -loc_fft_wmin*(tau_k - loc_fft_tmin)))).real() - 0.5*tail_coeff[1][s];
                    }
                }
            }
        }

        void FFT::transformTtoM(const ImTG& from, MatG& to, bool symmetric)
        {
            //bool symmetric = true;
            //fft_dt   = _beta/(tBins+2);
            //fft_tmin =  0.5*_beta/tBins;
            
            //TODO: this is a stupid workaround because the plan creation is not thread safe
            //TODO: for symmetric storage 2*nMF == tBins should be enforced
            const int tBins = static_cast<RealT>(_CONFIG_maxTBins);
            const int nMF = from.rows();
            const RealT tail_coeff1 = 1.;
            const RealT f1 = _beta/static_cast<RealT>(tBins);
            const int minMF = 0;//-_CONFIG_maxMatsFreq;
            fftw_complex fftw3_output[tBins];
            fftw_complex fftw3_input[tBins];
            auto planTtoM = fftw_plan_dft_1d(tBins,fftw3_input,fftw3_output,FFTW_BACKWARD,FFTW_ESTIMATE);
            const RealT loc_fft_dt = _beta/(tBins);
            const RealT loc_fft_tmin = 0;//0.5*loc_fft_dt;
            const RealT loc_fft_wmin = mFreqS( minMF, _beta);
            for(int s=0;s<2;s++){
                memset(fftw3_input, 0, tBins*sizeof(fftw_complex));
                for(int n=0;n<tBins;n++){
                    const RealT tau_k = loc_fft_dt*n + loc_fft_tmin;
                    ComplexT input_tmp = f1*std::exp( ComplexT(0.0,loc_fft_wmin * tau_k))*(from(n,s) + 0.5);
                    // ( t1(tail_coeff1,0.,tau_k) );
                    fftw3_input[n][0] = input_tmp.real();
                    fftw3_input[n][1] = input_tmp.imag();
                }
                //fftw3_input[tBins-1][0] = 0.;
                //fftw3_input[tBins-1][1] = 0.;
                fftw_execute(planTtoM);

                //RealT fac = -0.5*_beta*(from(0,s)+1.0)*from(nMF-1,s)/static_cast<RealT>(tBins) ;
                for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                    RealT mf =  mFreqS(n,_beta);
                    ComplexT tmp = (ComplexT(fftw3_output[n][0],fftw3_output[n][1]));
                    to(n,s) = (tmp +1.0/ComplexT(0.,mf));
                    //*std::exp(ComplexT(0.0,loc_fft_tmin*(mf-loc_fft_wmin)));
                }
            }
            //LOG(ERROR) << to;
            //exit(0);
        }


        void FFT::transformTtoM_test(const ImTG& from, MatG& to, bool symmetric)
        {
            //fft_dt   = _beta/(tBins+2);
            //fft_tmin =  0.5*_beta/tBins;
            
            //TODO: this is a stupid workaround because the plan creation is not thread safe
            //TODO: for symmetric storage 2*nMF == tBins should be enforced
            const int tBins = (symmetric+1)*static_cast<RealT>(_CONFIG_maxTBins);
            const int nMF = from.rows();
            fftw_complex fftw3_output[tBins];
            fftw_complex fftw3_input[tBins];
            auto planTtoM = fftw_plan_dft_1d(tBins,fftw3_input,fftw3_output,FFTW_BACKWARD,FFTW_ESTIMATE);
            const RealT loc_fft_dt = _beta/(tBins);
            const RealT loc_fft_tmin = 0.;//0.5*loc_fft_dt;
            const int n_min = 0 ;//symmetric ? 0 : -nMF/2.;
            const int n_max = n_min + nMF;
            const RealT loc_fft_wmin = mFreqS( n_min, _beta);
            for(int s=0;s<2;s++){
                memset(fftw3_input, 0, tBins*sizeof(fftw_complex));
                Eigen::VectorXd from2(tBins);
                if(symmetric)
                {
                    from2(0) = -0.5;
                    from2(tBins-1) = -0.5;
                    for(int t=1; t < _CONFIG_maxTBins; t++)
                    {
                        from2(2*t) = from(t,s);
                        from2(2*t-1) = from(t-1,s) + (from(t,s)-from(t-1,s))/2.;
                    }
                }
                //LOG(ERROR) << from2;
                //exit(0);
                else
                {
                    from2 = from.col(s);
                }
                for(int n=0;n<tBins;n++){
                    const RealT tau_k = loc_fft_dt*n + loc_fft_tmin;
                    ComplexT input_tmp = (std::exp( ComplexT(0.0,loc_fft_wmin * tau_k))) * from2(n);
                    fftw3_input[n][0] = input_tmp.real();
                    fftw3_input[n][1] = input_tmp.imag();
                }
                fftw_execute(planTtoM);

                for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                    RealT mf =  mFreqS(n + n_min,_beta);
                    ComplexT tmp = _beta*ComplexT(fftw3_output[n][0],fftw3_output[n][1])/static_cast<RealT>(tBins);
                    if(symmetric)
                    {
                        to(n,s) = tmp*std::exp(ComplexT(0.0,loc_fft_tmin*(mf-loc_fft_wmin)));
                    }
                    else
                    {
                        to((n+nMF/2)%nMF,s) = tmp*std::exp(ComplexT(0.0,loc_fft_tmin*(mf-loc_fft_wmin)));
                    }
                    //LOG(WARNING) << 1./mFreqS(2*n, _beta)  << " : " << to(n,s); 
                }
            }
        }

        void FFT::transformMtoT_naive(const MatG &from, ImTG &to) const
        {
            LOG(WARNING) << "untested";
            LOG(WARNING) << "this transofrmation method has not been updated for symmetric storage!"; 
            for(int t=0;t<_CONFIG_maxTBins;t++)
            {
                RealT tau = t*_beta/_CONFIG_maxTBins;
                for(int s=0;s<_CONFIG_spins;s++)
                {
                    ComplexT sum = 0.0;
                    for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                        const ComplexT tmp(std::exp(ComplexT(0.0,-mFreq(n,_beta)*tau)));
                        sum += tmp*(from(n,s) - 1.0);

                    }
                    to(t,s) = std::real(sum/_beta) - 0.5;
                }
            }
        }

        void FFT::transformTtoM_naive(const ImTG &from, MatG &to)//, bool symmetric)
        {
            const int tBins = _CONFIG_maxTBins*(1);
            const RealT dt = _beta/(tBins);
            const int n_min = 0;//-_CONFIG_maxMatsFreq-1;
            const int n_max = n_min + _CONFIG_maxTBins;//_CONFIG_maxMatsFreq;
            const int nMF = _CONFIG_maxMatsFreq;
            const int nh = _CONFIG_maxMatsFreq/2;
            for(int s=0;s<_CONFIG_spins;s++)
            {
                for(int n= n_min; n < n_max ; n++)
                {
                    RealT mf = mFreqS(n-_CONFIG_maxMatsFreq/2, _beta);
                    ComplexT sum(0.,0.);
                    for(int t=0;t<tBins;t++)
                    {
                        RealT tau = t*dt + 0.*dt/2.;
                        sum += std::exp(ComplexT(0.0,mf*tau))*(from(t)+0.5);

                    }
                    //(n+nh)%nMF
                    LOG(ERROR) << dt*sum - ComplexT(0., 1./mf);
                    //to(n-n_min, s) = ComplexT(0., std::imag(dt*sum));
                }
                exit(0);
            }
        }

        /*void FFT::transformTtoM_naive(GreensFct* gf) const
        {
            const RealT dtau = _beta/(_CONFIG_maxTBins-1);
            for(int s=0;s<_CONFIG_spins;s++)
            {
                for(int n=0;n<_CONFIG_maxMatsFreq;n++)
                {
                    int ng = n - (!gf->isSymmetric())*_CONFIG_maxMatsFreq/2.;
                    RealT wn = gf->getMFreq(ng);
                    ComplexT sum = 0.0;
                    for(int t=0;t<_CONFIG_maxTBins;t++)
                    {
                        RealT tau0 = (t+1)*_beta/(_CONFIG_maxTBins+2);
                        RealT tau1 = (t+2)*_beta/(_CONFIG_maxTBins+2);
                        ComplexT f0 = std::exp(ComplexT(0.0,wn*tau0))*(*gf)(tau0,s);
                        ComplexT f1 = std::exp(ComplexT(0.0,wn*tau1))*(*gf)(tau1,s);
                        sum += (f0+f1)/2.;

                    }
                    gf->setByMFreq(ng,s,dtau*sum);
                }
            }
        }*/
        
        /*void FFT::transformMtoT_naive(GreensFct* gf) const
        { 
            LOG(WARNING) << "debug function, use transformMtoT instead";
            if( std::abs(gf->getTailCoef(1, 0) - 1.)  > 0.01 || std::abs(gf->getTailCoef(1, 1)-1) >0.01 ) LOG(WARNING) << "Unexpected tail coefficient: " << gf->getTailCoef(1, 0) << ", "<<gf->getTailCoef(1, 1);
            ImTG new_gf(_CONFIG_maxTBins, _CONFIG_spins);
            const RealT tau_min = 0.;
            const RealT w_min = 0.;
            for(int s=0;s<_CONFIG_spins;s++)
            {
                for(int t=0;t<_CONFIG_maxTBins;t++)
                {
                    RealT tau = (t)*_beta/(_CONFIG_maxTBins);
                    ComplexT sum = ComplexT(0., 0.);
                    for(int n=0;n<_CONFIG_maxMatsFreq;n++)
                    {
                        RealT wn = mFreqS(n,_beta);
                        //ComplexT c = _beta*std::exp(ComplexT(0., tau_min*(wn - w_min)))/_CONFIG_maxTBins;
                    }
                    new_gf(t,s) = std::real(sum);
                }
            }
            gf->setByT(new_gf);
        }*/


        void FFT::conv(const ImTG& kernel, const ImTG& g, ImTG& res)
        {
            const int padded_size = kernel.rows() + g.rows() - 1;
            std::vector<RealT> inputFFT(padded_size);
            std::vector<ComplexT> outputFFT(padded_size/2 + 1);
            std::vector<ComplexT> convIn(padded_size/2 + 1);
            std::vector<RealT> convOut(padded_size);

            fftw_plan planR2C = fftw_plan_dft_r2c_1d(
                    padded_size,
                    reinterpret_cast<double*>(&inputFFT[0]),
                    reinterpret_cast<fftw_complex*>(&outputFFT[0]),
                    FFTW_ESTIMATE);
            fftw_plan planC2R = fftw_plan_dft_c2r_1d(
                    padded_size,
                    reinterpret_cast<fftw_complex*>(&convIn[0]),
                    reinterpret_cast<double*>(&convOut[0]),
                    FFTW_ESTIMATE);


            // TODO: use eigen::map to increase speed by avoiding copy operation
            Gnuplot gp1("gnuplot -persist");
            Gnuplot gp2("gnuplot -persist");
            Gnuplot gp3("gnuplot -persist");
            std::vector<std::pair<RealT,RealT>> resGP;

            for(int s=0;s<2;s++)
            {
                std::fill(inputFFT.begin(), inputFFT.end(), 0.0);
                for(int n = 0; n < kernel.rows(); n++)
                {
                    inputFFT[n] = kernel(n,s),0.0;
                }
                fftw_execute(planR2C);

                std::copy(outputFFT.begin(), outputFFT.end(), convIn.begin());

                std::fill(inputFFT.begin(), inputFFT.end(), 0.0);
                for(int n = 0; n < g.rows(); n++)
                {
                    inputFFT[n] = g(n,s);
                }
                fftw_execute(planR2C);

                for(int n = 0; n < convIn.size(); n++)
                {
                    convIn[n] = convIn[n] * outputFFT[n];
                }
                fftw_execute(planC2R);

                for(int n = 0; n < res.rows(); n++)
                {
                    resGP.push_back( std::make_pair(n,convOut[n]/padded_size) );
                    res(n,s) = g(n,s) - (convOut[2*n]+convOut[2*n+1])/padded_size;
                }
                if(s==0)
                {
                    gp1 << "set output ' res.png'\n";
                    gp1 << "plot" << gp1.file1d( resGP ) << "with points title 'FFT convolution'" << std::endl;

                }
            }
            fftw_destroy_plan(planR2C);
            fftw_destroy_plan(planC2R);
        }

        /*  for(int s=0;s<2;s++){
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
            }*/

        std::vector<RealT> FFT::conv_naive(const ImTG& g0, const ImTG& kernel, const RealT totalSign)
        {
            std::vector<RealT> res(2*g0.rows());
            //LOG(INFO) << "g0 rows, kernel rows: " << g0.rows() << ", " << kernel.rows();
            //std::cout << g0 << std::endl << "g0---conv---kernel\n" << kernel<< std::endl;
            for(int s = 0; s < 2; s++)
            {

                for(int i = 0; i < g0.rows(); i++)
                {
                    RealT tau = _beta*static_cast<RealT>(i)/g0.rows();
                    RealT tmp = 0.0;
                    for(int j = 0; j < kernel.rows(); j++)
                    {
                        if(kernel(j,s) == 0.0) continue;
                        RealT taup = _beta*static_cast<RealT>(j)/kernel.rows(); 
                        tmp += static_cast<int>(2*(tau>=taup)-1)*g0(static_cast<int>((tau-taup + static_cast<int>(tau<taup)*_beta)*g0.rows()/_beta),s)*kernel(j,s);
                        //LOG(INFO) << "tau, taup \t" <<  tau << ", " << taup << "  \t tmp: " << tmp;
                    }
                    //LOG(INFO) << "tau_conv: " << tau << ", sum: " << tmp << ", tmp/totalSign: " <<  tmp/totalSign << ", g0: " << g0(i,s) << ". setting: " << g0(i,s) - tmp/totalSign;
                    res[g0.rows()*(s == UP) + i] = g0(i,s) - tmp/totalSign;

                }

            }
            return res;
        }

        /*
           void _TEST_naiveImpl(fftw_complex* input,fftw_complex* output, RealT _beta)
           {

           for(int i=0;i<_CONFIG_maxTBins;i++){
           for(int s=0;s<_CONFIG_spins;s++){
           output[arrIndex(i,s)][0] = 0.0;
           output[arrIndex(i,s)][1] = 0.0;
           const RealT tau = (i/_CONFIG_maxTBins);
           for(int n=0;n<_CONFIG_maxMatsFreq;n++){
           const ComplexT gwn = ComplexT(input[arrIndex(n,s)][0],input[arrIndex(n,s)][1]);

           output[arrIndex(i,s)][0] += (gwn*std::exp( ComplexT(0.0,-mFreq(n,_beta)* tau) )).real();
           output[arrIndex(i,s)][1] += (gwn*std::exp( ComplexT(0.0,-mFreq(n,_beta)* tau) )).imag();
           }
           output[arrIndex(i,s)][0] = output[arrIndex(i,s)][0]/_beta;
           output[arrIndex(i,s)][1] = output[arrIndex(i,s)][1]/_beta;
           }
           }
           }

           void _TEST_naive2Impl(fftw_complex* input,fftw_complex* output, RealT _beta)
           {

           for(int i=0;i<_CONFIG_maxTBins;i++){
           for(int s=0;s<_CONFIG_spins;s++){
           output[arrIndex(i,s)][0] = 0.0;
           output[arrIndex(i,s)][1] = 0.0;
           const RealT tau = (i/_CONFIG_maxTBins);

           for(int n=0;n<_CONFIG_maxMatsFreq;n++){
           const ComplexT gwn = ComplexT(input[arrIndex(n,s)][0],input[arrIndex(n,s)][1]);

           output[arrIndex(i,s)][0] += ( gwn*std::exp(ComplexT(0.0,-mFreq(n,_beta)*2.0* M_PIl*tau/_beta)) ).real();
           output[arrIndex(i,s)][1] += ( gwn*std::exp(ComplexT(0.0,-mFreq(n,_beta)*2.0* M_PIl*tau/_beta)) ).imag();
           }
           ComplexT tmp(output[arrIndex(i,s)][0], output[arrIndex(i,s)][1]);
           output[arrIndex(i,s)][0] = (std::exp(ComplexT(0.0,-tau*M_PIl/_beta))*tmp/_beta).real();
           output[arrIndex(i,s)][1] = (std::exp(ComplexT(0.0,-tau*M_PIl/_beta))*tmp/_beta).imag();
           }
           }
           }*/

} //end namespace DMFT

