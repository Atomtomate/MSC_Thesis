#ifndef FFT_HPP_
#define FFT_HPP_

#include "Config.hpp"


#include "fftw3.h"
#include <algorithm>

//TODO: use fftw3-mpi and/or GPU fft
//#include "gnuplot-iostream.h"

namespace DMFT
{
    class FFT
    {

        private:
            //TODO: check for _CONFIG_maxMatsFreq != _CONFIG_maxTBins
            RealT _beta;
            fftw_complex fftw3_input[_CONFIG_maxMatsFreq];
            fftw_complex fftw3_output[_CONFIG_maxMatsFreq];


            fftw_plan planMtoT;
            fftw_plan planTtoM;
            RealT fft_tmin;
            RealT fft_dw;
            RealT fft_dt;
            RealT fft_wmin;

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

            /*! performs a convolution of f and g using fft
             *  fp = FFT(f), gp = FFT(g) -> iFFT(fp * gp): the result is int f(t' - t) g(t) dt 
             *
             *  @param  [in]  kernel   function 1
             *  @param  [in]  g   function 2
             *
             *  @param  [out] res g(t') - int g(t' - t) kernel(t) dt
             */
            void conv(const ImTG& kernel, const ImTG& g, ImTG& res);

            /*! performs naive convolution
             *
             *  @param [in] g0
             *  @param [in] s  kernel
             *
             *  @return     g0(t) - integral g0(t')s(t-t') dt' as vector, res_DOWN[...], res_UP[...]
             */
            std::vector<RealT> conv_naive(const ImTG& g0, const ImTG& kernel, const RealT totalSign);

            /*! FFT from imaginary time to Matsubara space.
             *  TODO: TRIQS ref, detailed description
             *
             *  @param [in] from    ImTGF as Eigen::ArrayXXd
             *  @param [out] to     MatGF as Eigen::ArrayXXcd
             */
            void transformTtoM(const ImTG& from, MatG& to);

    };

}   //end namespace DMFT
#endif
