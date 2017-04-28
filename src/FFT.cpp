#include "FFT.hpp"
namespace DMFT
{

    /*  void FFT::fft(const MatG &from, MatG &to)
        {

        }

        void FFT::inv_fft(const MatG &from, MatG &to)
        {

        }
        */

    //TODO: use http://www.fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html
    void FFT::transformMtoT(const MatG &from, ImTG &to)
    {
        const RealT absErr = 0.01;
        const int shift = static_cast<int>((_CONFIG_maxTBins - _CONFIG_maxMatsFreq)/2);
        for(int s=0;s<_CONFIG_spins;s++){
            memset(fftw3_input, 0, _CONFIG_maxTBins*sizeof(fftw_complex));
            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                const RealT wn_k = fft_wmin + n*fft_dw;
                ComplexT input_tmp = std::exp( ComplexT(0.0,-fft_tmin * wn_k ))*( from(n,s) - 1.0/(ComplexT(0.0,mFreq(n,_beta))) )/_beta;
                fftw3_input[n][0] = input_tmp.real();
                fftw3_input[n][1] = input_tmp.imag();
            }
            fftw_execute(planMtoT);

            for(int n=0;n<_CONFIG_maxTBins;n++){
                const RealT tau_k = fft_dt*n;
                ComplexT tmp = ComplexT(fftw3_output[n][0],fftw3_output[n][1]);
                tmp = tmp*std::exp(ComplexT(0.0,-fft_wmin*tau_k));
                fftw3_output[n][0] = tmp.real() - 0.5; // -0.5 FT of - 1.0/(ComplexT(0.0,mFreq(n,_beta)))
                fftw3_output[n][1] = tmp.imag();
                //if(fftw3_output[n][1] > absErr) LOG(WARNING) << "imaginary part of im T Greens function has significant contribution";
                //TODO: different sizs of t and n bins??
                to(n,s) = fftw3_output[n][0];
            }
            //TODO: test std::copy
        }
    }

    void FFT::transformTtoM(const ImTG& from, MatG& to)
    {
        for(int s=0;s<_CONFIG_spins;s++){
            memset(fftw3_input, 0, _CONFIG_maxTBins*sizeof(fftw_complex));
            for(int n=0;n<_CONFIG_maxTBins;n++){
                const RealT tau_k = fft_dt*n + fft_tmin;
                ComplexT input_tmp = (std::exp( ComplexT(0.0,fft_wmin * tau_k))) * from(n,s);
                fftw3_input[n][0] = input_tmp.real();
                fftw3_input[n][1] = input_tmp.imag();
            }
            fftw_execute(planTtoM);

            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                const RealT mf =  mFreq(n,_beta);
                ComplexT tmp = ComplexT(fftw3_output[n][0],fftw3_output[n][1])*ComplexT(_beta/_CONFIG_maxTBins,0.0);
                tmp = tmp*std::exp(ComplexT(0.0,-fft_tmin*mf));
                fftw3_output[n][0] = tmp.real();
                fftw3_output[n][1] = tmp.imag();
                to(n,s) = tmp;
            }
        }
    }

    void FFT::transformMtoT_naive(const MatG &from, ImTG &to)
    {
        for(int t=0;t<_CONFIG_maxTBins;t++)
        {
            RealT tau = t*_beta/_CONFIG_maxTBins;
            for(int s=0;s<_CONFIG_spins;s++)
            {
                ComplexT sum = 0.0;
                for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                    const ComplexT tmp(std::exp(ComplexT(0.0,-mFreq(n,_beta)*tau)));
                    sum += tmp*(from(n,s) - 1.0/ComplexT(0.0,mFreq(n,_beta)));

                }
                to(t,s) = std::real(sum/_beta) - 0.5;
            }
        }
    }

    void FFT::transformTtoM_naive(const ImTG &from, MatG &to)
    {
        //TODO: implement
        for(int t=0;t<_CONFIG_maxTBins;t++)
        {
            RealT tau = t*_beta/_CONFIG_maxTBins;
            for(int s=0;s<_CONFIG_spins;s++)
            {
                ComplexT sum = 0.0;
                for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                    const ComplexT tmp(std::exp(ComplexT(0.0,-mFreq(n,_beta)*tau)));
                    sum += tmp*(from(n,s) - 1.0/ComplexT(0.0,mFreq(n,_beta)));

                }
                to(t,s) = std::real(sum/_beta) - 0.5;
            }
        }
    }


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

