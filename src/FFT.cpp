#include "FFT.hpp"
namespace DMFT
{


void FFT::transformMtoT(const MatG &from, ImTG &to)
{
    const RealT absErr = 0.01;
    for(int s=0;s<_CONFIG_spins;s++){
	for(int n=0;n<_CONFIG_maxMatsFreq;n++){
	    const RealT wn_k = fft_wmin + n*fft_dw;
	    ComplexT input_tmp = std::exp( ComplexT(0.0,-fft_tmin * wn_k ))*( from(n,s) - 1.0/(ComplexT(0.0,mFreq(n,_beta))) )/_beta;
	    fftw3_input[n][0] = input_tmp.real();
	    fftw3_input[n][1] = input_tmp.imag();
	}
	fftw_execute(planMtoT);

	for(int n=0;n<_CONFIG_maxMatsFreq;n++){
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

void FFT::transformTtoM(const ImTG &from, MatG &to)
{
    for(int s=0;s<_CONFIG_spins;s++){
	for(int n=0;n<_CONFIG_maxTBins;n++){
	    const RealT tau_k = fft_dt*n + fft_tmin;
	    ComplexT input_tmp = (std::exp( ComplexT(0.0,fft_wmin * tau_k))) * from(n,s);
	    fftw3_input[n][0] = input_tmp.real();
	    fftw3_input[n][1] = input_tmp.imag();
	}
    }
    fftw_execute(planTtoM);

    for(int s=0;s<_CONFIG_spins;s++){
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
	for(int t=0;t<_CONFIG_maxTBins;t++){
        RealT tau = t*_beta/_CONFIG_maxTBins;
        for(int s=0;s<_CONFIG_spins;s++){
            ComplexT sum = 0.0;
	        for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                const ComplexT tmp(std::exp(ComplexT(0.0,-mFreq(n,_beta)*tau)));
                sum += tmp*(from(n,s) - 1.0/ComplexT(0.0,mFreq(n,_beta)));

            }
            to(t,s) = std::real(sum/_beta) - 0.5;
        }
    }
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

