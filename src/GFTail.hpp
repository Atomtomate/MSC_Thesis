#ifndef GFTAIL_HPP_
#define GFTAIL_HPP_

#include "Config.hpp"
#include "Constants.h"
#include "FFT.hpp"
//#include "GreensFct.hpp"

#include <vector>
#include <complex>
#include <iostream>
#include <string>

#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>

namespace DMFT
{

    struct GFTail
    {
        // no function pointer because we want to accept lambdas which capture
        // for example beta
        //(int)(_CONFIG_maxMatsFreq/30)
        GFTail() : nC(5), first( std::max(_CONFIG_maxMatsFreq/20, 20) ), last(_CONFIG_maxMatsFreq-10), fitFct([](int n, int i, RealT beta) { return std::pow(1.0/DMFT::ComplexT(0.,mFreqS(n, beta)),i);}) {}
        std::function<ComplexT(int, int, RealT)> fitFct;                               
        unsigned int first;
        unsigned int last;
        unsigned int nC;
    };

    // return std::pow(1.0/DMFT::mFreqS(n, beta),i);
    //DMFT::GFTail<8, 1000, DMFT::ComplexT> gft(g0, [&beta](int n, int i) { return std::pow(1.0/DMFT::ComplexT(0.,mFreqS(n, beta)),i);}); //
    //DMFT::GFTail<8, 1000, DMFT::ComplexT> gft(g0, 10. ); //
    //gft.fitTail(24);

    //template<unsigned char nC, unsigned int nSamples, typename T, typename F, typename TailType = typename std::result_of<F&(T,int)>::type>
    /*template< unsigned char nC, unsigned int nSamples, typename T>
    class GFTail
    {
        public:
            //GFTail(const GreensFct &gf, std::function<T(int, int)> f):
            //    gf(gf), config(config), func(f)
            //{
            //     A.resize(nSamples, nC);
                 //coeffs(nC);
                 //rhs(nSamples);
            //}

            GFTail(const GreensFct &gf, const RealT beta):
                gf(gf), config(config), func([&beta](int n, int i) { return std::pow(1.0/DMFT::ComplexT(0.,mFreqS(n, beta)),i);})
            {
                 A.resize(nSamples, nC);
            }


            void fitTail(const unsigned int usableFrom)
            {
                for(unsigned char f = 1; f < _CONFIG_spins; f++)
                {
                    for(int i = 0; i < nSamples && (i + usableFrom < _CONFIG_maxMatsFreq); i++)
                    {
                        if(i+usableFrom > _CONFIG_maxMatsFreq) LOG(ERROR) << "too many values for tail fit requested";
                        for(int j = 0; j < nC; j++)
                        {
                            //LOG(INFO) << usableFrom+i << " : "<< gf.getByMFreq(usableFrom+i, f) << " : "<<  gf.getMFreq(usableFrom + i) << " : " << func( gf.getMFreq(usableFrom + i) ,j);
                            //A(i,j) = func( usableFrom + i ,j);
                            //rhs(i) = std::imag(gf.getByMFreq(usableFrom+i, f));
                            A(i,j) = func( usableFrom + i ,j);
                            rhs(i) = gf.getByMFreq(usableFrom+i, f), T();
                        }
                    }
                    coeffs = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs).real().transpose();
                    if(gf.isSymmetric())
                    {
                        for(int j = 0; j < nC; j+=2)
                            coeffs[j] = 0.;
                    }
                }
            }



        private:
            //Eigen::Matrix<T, nSamples, nC>
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
            Eigen::Matrix<T, nSamples, 1> rhs;
            //Eigen::Matrix<T, Eigen::Dynamic, 1> rhs;
            Eigen::Matrix<RealT, 1, nC> coeffs;
            //Eigen::Matrix<T, 1, Eigen::Dynamic> coeffs;
            const std::function<T(int, int)> func;
            const GreensFct &gf; 
            const Config &config; 
    };
    */

/*
        */
}

#endif
