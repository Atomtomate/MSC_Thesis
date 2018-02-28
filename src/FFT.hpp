#ifndef FFT_HPP_
#define FFT_HPP_

#include "Config.hpp"

#include "fftw3.h"
#include <algorithm>

//TODO: use fftw3-mpi and/or GPU fft
//#include "gnuplot-iostream.h"
namespace DMFT
{
    class GreensFct;
    class FFT
    {

        private:
            //TODO: check for _CONFIG_maxMatsFreq != _CONFIG_maxTBins
            RealT _beta;
            //TODO: use r2c and c2r
            fftw_complex fftw3_input[_CONFIG_maxTBins];
            fftw_complex fftw3_output[_CONFIG_maxTBins];

            fftw_plan planMtoT;
            fftw_plan planTtoM;
            RealT fft_tmin;
            RealT fft_dw;
            RealT fft_dt;
            RealT fft_wmin;

        //TODO: this can all be obtained from GrennFct object directly
            ComplexT tail(const RealT wn, const std::vector<std::array<RealT,2> >& tail, const int spin) const;

        //TODO: this can all be obtained from GrennFct object directly
            RealT ftTail(const RealT tau, const std::vector<std::array<RealT,2> >& tail, const int spin) const;

        public:
            FFT(const RealT beta): _beta(beta)
        {
            planMtoT = fftw_plan_dft_1d(_CONFIG_maxTBins,fftw3_input,fftw3_output,FFTW_FORWARD,FFTW_ESTIMATE);
            planTtoM = fftw_plan_dft_1d(_CONFIG_maxTBins,fftw3_input,fftw3_output,FFTW_BACKWARD,FFTW_ESTIMATE);
            fft_wmin = mFreq(0,_beta);
            fft_dt   = _beta/static_cast<RealT>(_CONFIG_maxTBins);      // for FFT matsFreq == TBins
            fft_dw   = 2.0*boost::math::constants::pi<RealT>()/_beta;
            fft_tmin = 0;//fft_dt/2.0;
        }

            virtual ~FFT()
            {
                fftw_destroy_plan(planMtoT);
                fftw_destroy_plan(planTtoM);
            }

            /*! FFT from Matsubara space t>0 imaginary time.
             *  TODO: TRIQS ref, detailed description
             *
             *  @param [in] to      MatGF as Eigen::ArrayXXcd
             *  @param [out] from   ImTGF as Eigen::ArrayXXd
             */
            void transformMtoT(const MatG &from, ImTG &to, bool symmetric);

            /*! FFT from Matsubara space t>0 imaginary time.
             *  TODO: TRIQS ref, detailed description
             *
             *  @param [in] to      MatGF as Eigen::ArrayXXcd
             *  @param [in] tail    std::vector<std::array<RealT,2> > with analytic GF tail coefficients
             *  @param [out] from   ImTGF as Eigen::ArrayXXd
             */
            void transformMtoT(const MatG& from, ImTG& to, const std::vector<std::array<RealT,2> >& tail, bool symmetric);

            void transformMtoT_naive(const MatG& from, ImTG& to) const;
            void transformMtoT_naive(GreensFct* gf) const;

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

            void transformTtoM_naive(const ImTG &from, MatG &to) const;
            void transformTtoM_naive(GreensFct* gf) const;
    };

}   //end namespace DMFT
#endif
