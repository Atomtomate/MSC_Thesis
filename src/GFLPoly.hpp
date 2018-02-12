#ifndef GF_L_POLY_HPP_
#define GF_L_POLY_HPP_


#include "Config.hpp"
#include "GreensFct.hpp"

#include <initializer_list>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>

namespace DMFT
{

    class GFLPoly
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            GFLPoly(GreensFct* const gf, const Config& conf): gf(gf), conf(conf)
            {
                ComplexT mo(-1.,0.);
                ComplexT i(0.,1.);
                Tnl.resize(_CONFIG_maxMatsFreq, _CONFIG_maxLPoly);
                for(int n=0; n<_CONFIG_maxMatsFreq; n++)
                {
                    RealT x(0.5*conf.beta*mFreqS(n, conf.beta));
                    for(int l = 0; l< _CONFIG_maxLPoly; l++)
                    {
                        Tnl(n,l) = std::pow(mo, n)*std::pow(i,l+1)*std::sqrt(2*l+1)*boost::math::sph_bessel(l, x);
                    }
                }
            }

            inline void setGl(const unsigned int l, const int spin, const RealT val)
            {
                Gl[spin](l) = val;
            }

            inline ComplexT getByM(const int n, const int spin) const
            {
                return (Tnl.row(n)).dot(Gl[spin]);
            }

            MatG getByM(const int spin) const
            {
                return Tnl*Gl[spin];
            }

            RealT getByT(RealT tau, int spin)
            {
                RealT res = 0;
                RealT x = 2*tau/conf.beta - 1;
                for(int l=0; l < _CONFIG_maxLPoly; l++)
                {
                    res += boost::math::legendre_p(l,x)*Gl[spin](l)*std::sqrt(2*l+1);
                }
                return res/conf.beta;
            }

            void computeTail(void)
            {
                LOG(ERROR) << "tail computation not yet implemented";
                std::array<RealT,_CONFIG_tailOrder> coeffs{};
                for(int l=0; l<_CONFIG_maxLPoly;l++)
                {
                    if(l%2 == 0)
                    {
                        coeffs[1] = - 2.*std::sqrt(2.*l + 1)*Gl[0](l)/conf.beta;
                        coeffs[3] = - (l+2.)*(l+1.)*(l-1.)*std::sqrt(2.*l + 1)*Gl[0](l)/(conf.beta*conf.beta*conf.beta);
                    }
                    else
                    {
                        coeffs[2] = 2.*l*(l+1.)*std::sqrt(2.*l + 1)*Gl[0](l)/(conf.beta*conf.beta);
                    }
                }
                LOG(ERROR) << coeffs[1] << ", " << coeffs[2] << ", " << coeffs[3];
            }

            void setGF(void)
            {
                ImTG g_it = ImTG::Zero(_CONFIG_maxTBins, _CONFIG_spins);
                MatG g_mf = MatG::Zero(_CONFIG_maxMatsFreq, _CONFIG_spins);
                for(int f = 0; f<_CONFIG_spins; f++)
                {
                    //g_mf.col(f) = getByM(f);
                    for(int n=0; n<_CONFIG_maxMatsFreq; n++)
                    {
                        g_mf(n,f) = Tnl(n,0)*Gl[f](0);
                        for(int l = 1; l<_CONFIG_tailOrder;l++)
                            g_mf(n, f) += Tnl(n,l)*Gl[f](l);
                    }
                    for(int t=0; t<_CONFIG_maxTBins;t++)
                    {
                       g_it(t,f) = getByT(conf.beta*t/_CONFIG_maxTBins, f);
                    }
                }
                computeTail();
                gf->setByMFreq(g_mf);
                gf->setByT(g_it);
                gf->markMSet();
                gf->markTSet();
            }


        private:
            GreensFct * gf;
            const Config& conf;
            std::array<Eigen::Matrix<RealT, _CONFIG_maxLPoly,1>,_CONFIG_spins> Gl;
            CMatrixT Tnl;
    };
}
#endif
