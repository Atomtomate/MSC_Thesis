#ifndef FFT_HPP_
#define FFT_HPP_

#include "Config.hpp"
#include "fftw3.h"

//#include "gnuplot-iostream.h"

namespace DMFT
{
class FFT
{

public:
    FFT(const RealT beta): _beta(beta)
    {
	    planMtoT = fftw_plan_dft_1d(_CONFIG_maxMatsFreq,fftw3_input,fftw3_output,FFTW_FORWARD,FFTW_ESTIMATE);
	    planTtoM = fftw_plan_dft_1d(_CONFIG_maxMatsFreq,fftw3_input,fftw3_output,FFTW_BACKWARD,FFTW_ESTIMATE);
	    fft_wmin = mFreq(0,_beta);
	    fft_dt   = _beta/static_cast<RealT>(_CONFIG_maxTBins);      // for FFT matsFreq == TBins
	    fft_dw   = 2.0*boost::math::constants::pi<RealT>()/_beta;
	    fft_tmin = fft_dt/2.0;
    }

    virtual ~FFT()
    {
	    fftw_destroy_plan(planMtoT);
	    fftw_destroy_plan(planTtoM);
    }

    /*! FFT from Matsubara space t>o imaginary time.
     *  TODO: TRIQS ref, detailed description
     *
     *  @param [in] to      MatGF as Eigen::ArrayXXcd
     *  @param [out] from   ImTGF as Eigen::ArrayXXd
     */
    void transformMtoT(const MatG& from, ImTG& to);

    void transformMtoT_naive(const MatG& from, ImTG& to);

    /*! FFT from imaginary time to Matsubara space.
     *  TODO: TRIQS ref, detailed description
     *
     *  @param [in] from    ImTGF as Eigen::ArrayXXd
     *  @param [out] to     MatGF as Eigen::ArrayXXcd
     */
    void transformTtoM(const ImTG& from, MatG& to);

private:
    RealT _beta;
	fftw_complex fftw3_input[_CONFIG_maxMatsFreq];
	fftw_complex fftw3_output[_CONFIG_maxMatsFreq];
	fftw_plan planMtoT;
	fftw_plan planTtoM;
	RealT fft_tmin;
	RealT fft_dw;
	RealT fft_dt;
	RealT fft_wmin;
};

}   //end namespace DMFT
#endif
