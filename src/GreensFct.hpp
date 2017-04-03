#ifndef GREENS_FCT_HPP_
#define GREENS_FCT_HPP_

#include "Config.hpp"
#include "Constants.h"
#include "FFT.hpp"

#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <sstream>

namespace DMFT
{


    /*! Implementation for functionality of the Green's funtion  
     *  
     *  \f[ \Delta(\sigma,w_n,\epsilon_k,V(k,\sigma )) = \sum_k \|V(k,\sigma )\|^2 /(i*\omega_n - \epsilon_k) \f]	(8.14)
     *  \f[ G0^{-1} = i*\omega_n + \mu - 0.5*U - \Delta(\omega_n) \f] 	                                        (8.15)
     *  \f[ [M^{-1}]_{ij} = G_0 - \alpha * \delta_ij   \f]                                             		(8.34)
     *
     *  Implementation details:
     *  - alignment \f$ G_\downarrow (\tau_1, \tau_2, ...) , G_\uparrow (\tau_1, \tau_2, ...) \f$
     *
     *  TODO: finish implementation of GreensFct as expression template
     *  TODO: remember to add site dependency for hubbard model
     */
    class GreensFct//: public VecExpression<GreensFct>
    {

        const int MAX_T_BINS = _CONFIG_maxTBins;
        const int MAX_M_FREQ = _CONFIG_maxMatsFreq;
        const int SPINS	= _CONFIG_spins;

        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        GreensFct(RealT beta): beta(beta), fft(beta)
        {
            g_it = ImTG::Zero(MAX_T_BINS, SPINS);
            g_wn = MatG::Zero(MAX_M_FREQ, SPINS);
        }

        GreensFct(GreensFct const&) = delete; // TODO: implement

        virtual ~GreensFct(){}

        /*! Call this after setting \f$ G(\tau) \f$ 
         *  @return 1 if previously set, 0 otherwise
         */
        int markTSet(void){return tSet++;}

        /*! Call this after setting \f$ G(i \omega_n) \f$.
         *  @return returns 1 if previously set, 0 otherwise
         */
        int markMSet(void){return mSet++;}

        /*! Sets G(i \omega_n) = val
         *
         *  @param [in] n Matsubara frequency 
         *  @param [in] spin \f$ \sigma \in N \f$
         *  @param [in] val value of \f$ G_\sigma (\tau ) \f$
         */
        inline void setByMFreq(const int n, const int spin, const ComplexT val)	{ g_wn(n,spin) = val; }

        /*! Sets G(i \omega_n) = val
         *
         *  @param [in] g Matsubara Green's function (as Eigen::ArrayXXd) 
         */
        void setByMFreq(MatG &g)	{ g_wn = g; mSet = 1;}

        /*! Sets \f$ G(t) = val, G(0) = val \Rightarrow G(0^+) = val \f$
         *
         *  @param [in] t imaginary time \f$ \tau \f$
         *  @param [in] spin \f $\sigma \in N \f$
         *  @param [in] val value of Green's function for \sigma at \tau
         */
        inline void setByT(const RealT t, const int spin, const RealT val){
            const unsigned int i = tIndex(t);
            //if(t==1) g_it(0,spin) = val;	// NOTE: redudant since tIndex shifts up from 0, remove?
            g_it(i,spin) = val;
        }
        void setByT(ImTG &g) {g_it = g; tSet = 1;}

        /*! get Matsubara Green's function \f$ G_\sigma(i \omega_n) \f$ by n.
         *  @param  n Matsubara frequency \f$ \omega_n \f$
         *  @param  spin \f$ \sigma \in N \f$
         *
         *  @return
         */
        inline ComplexT getByMFreq(const int n, const unsigned int spin) const		{ return g_wn(n,spin); }

        /*! @brief	imaginary time Green's function \f$ G_\sigma(\tau) \f$
         *  		storage order is t first (inner loop over t)
         *  
         *  @param [in] t imaginary time \f$ \tau \f$
         *  @param [in] spin \f$ \sigma \in N \f$
         *
         *  @return G_sigma(  \f$ \tau \f$ )
         */
        inline RealT operator() (const RealT t, int unsigned spin) const{
            const int sgn = 2*(t>0)-1;
            //if(t<0){ t+=beta; return -1.0*g_it(tIndex(t), spin);}
            return sgn*g_it( tIndex( t + (t<0)*beta ), spin);
        }
        template <class T> RealT operator()(T) = delete; // C++11

        /*! Return space separated string for the Matsubara Green's function.
         *  Format: \f$ \omega_n \quad Re[G_\downarrow (i \omega_n)] \quad  Im[G_\downarrow (i \omega_n)] \quad Re[G_\uparrow (i \omega_n)] \quad  Im[G_\uparrow (i \omega_n)]\f$ 
         *  
         *  @return string with data.
         */
        std::string getMGFstring(void) const
        {
            std::stringstream res;
            res << std::fixed << std::setw(8)<< "mFreq \tSpin0 Re\tSpin0 Im\tSpin1 Re\tSpin1 Im\n";
            for(int n=0; n < MAX_M_FREQ; n++)
            {
                res << std::setprecision(8) << mFreq(n, beta) << "\t" <<  g_wn(n,0).real() << "\t" \
                    << g_wn(n,0).imag() << "\t" << g_wn(n,1).real() << "\t" << g_wn(n,1).imag() << "\n" ;
            }
            return res.str();
        }

        std::string getMaxEntString(void) const
        {
            std::stringstream res;
            res << std::fixed << std::setw(8);
            for(int n=static_cast<int>(MAX_M_FREQ/2 ); n < MAX_M_FREQ; n++)
            {
                res << std::setprecision(8) << mFreq(n,beta) << " " << g_wn(n,0).imag() << 0.0001 << "\n" ;
            }
            return res.str();
        }

        /*! Return space separated string for the Matsubara Green's function.
         *  Format: \f$ \omega_n \quad Re[G_\downarrow (\tau)] \quad  Im[G_\downarrow (\tau )] \quad Re[G_\uparrow (\tau )] \quad  Im[G_\uparrow (\tau )]\f$ 
         *  
         *  @return string with data.
         */
        std::string getITGFstring(void) const
        {
            std::stringstream res;
            res << std::fixed << std::setw(8)<< "iTime \tSpin up \tSpin down\n";
            for(int i =0; i < MAX_T_BINS; i++)
            {
                res << std::setprecision(8) << std::setw(8) << (beta*i)/MAX_T_BINS << "\t"\
                    << g_it(i,0)<< "\t" << g_it(i,1) <<"\n";
            }
            return res.str();
        }

        /*! Does a FFT from Matsubara to imaginary time Green's function and stores it in internal storage.
        */
        inline void transformMtoT(void){ fft.transformMtoT_naive(g_wn, g_it); markTSet(); }

        /*! Does a FFT from Matsubara to imaginary time Green's function and stores it in target.
         *  @param	[out]	Eigen::ArrayXXd where imaginary time Green's function will be stored
         */
        inline void transformMtoT(ImTG& target){ fft.transformMtoT_naive(g_wn, target); markTSet(); }

        /*! Does a FFT from imaginary time to Matsubara Green's function and stores it.
        */
        inline void transformTtoM(void){ fft.transformTtoM(g_it, g_wn); markMSet(); }   

        /*! Does a FFT from imaginary time to Matsubara Green's function and stores it in target.
         *  @param	[out]	Eigen::ArrayXXcd where Matsubara Green's function will be stored
         */
        inline void transformTtoM(MatG& target){ fft.transformTtoM(g_it, target); markMSet(); }

        /*! Compute \f[ ()G(i \omega_n)^{-1} - s)^{-1}  \f] and also transform to \f[ G(\tau) \f].
         *  @param  [in]  s shift
         */
        inline void shift(const RealT s) { g_wn = ((g_wn.cwiseInverse() - s).cwiseInverse()).eval() ; markMSet(); transformMtoT();}

        const MatG& getMGF() const { return g_wn; }
        const ImTG& getItGF() const { return g_it; }

        // TODO: overload operators and use expression templates



        void setParaMagnetic(void)
        {
            /*LOG(DEBUG) << g_wn;
              LOG(DEBUG) << "===";
              LOG(DEBUG) << g_wn.leftCols(1);
              LOG(DEBUG) << "===";
              LOG(DEBUG) << g_wn.rightCols(1);
              LOG(DEBUG) << "---";*/
            g_wn.rightCols(1) = 0.5*(g_wn.leftCols(1) + g_wn.rightCols(1)).eval();
            g_wn.leftCols(1) = g_wn.rightCols(1).eval();
            transformMtoT();
        }
        protected:
        RealT beta;
        //TODO ASSERT MATSFREQ == TBINS !!!!!
        ImTG	g_it; // col major -> spin outer loop
        MatG	g_wn; // col major -> spin outer loop
        FFT 	fft;
        //Config& conf;

        int mSet = 0;	// track whether matsuabra or iTime GF has been set
        int tSet = 0;	// TODO implement

        /*! computes index from imaginary time. G(0) => G(0^+), G(beta) => G(beta^-)
         *  @param [in] t imaginary time \tau
         *
         *  @return index for lookup in g_it
         */
        inline int tIndex(const RealT t) const {
            const int res = static_cast<int>(t*static_cast<RealT>(MAX_T_BINS)/beta);
            return res+(res == 0);
        }


    };

    //using GreensFct = GreensFct;	// useless for now, may be useful for later changes to GreensFct
    //using GreensFct = GreensFct<_CONFIG_maxTBins, _CONFIG_maxMatsFreq, _CONFIG_spins>;
}	//end namespace DMFT

#endif
