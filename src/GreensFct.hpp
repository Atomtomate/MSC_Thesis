#ifndef GREENS_FCT_HPP_
#define GREENS_FCT_HPP_

#include "Config.hpp"
#include "Constants.h"
#include "FFT.hpp"
#include "GFTail.hpp"

#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <sstream>

namespace DMFT
{
//DMFT::ComplexT c_fit_fct(int n, int i)
//{
//    return std::pow(1./DMFT::ComplexT(0., DMFT::mFreqS(n, 10.)),i);
//}

    /*ComplexT defaultTail(int n, int i)
    {
        return std::pow(1.0/ComplexT(0.,mFreqS(n, 10.)),i);
    }*/
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
     *  TODO: use cubic spline interpolation for FFT - halfway done
     *  TODO: high freq. expansion (can use cubic spline)
     *  TODO: convert symmetric to template paramter when C++17 becomes standard
     *  TODO: remember to add site dependency for hubbard model
     */
    class GreensFct//: public VecExpression<GreensFct>
    {

        //ComplexT (*tailFunc)(int, int);
        static constexpr int MAX_T_BINS = _CONFIG_maxTBins;
        static constexpr int MAX_M_FREQ = _CONFIG_maxMatsFreq;
        static constexpr int SPINS	= 2;
        static constexpr int TAIL_ORDER = 7;

        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            GreensFct(RealT beta, const bool symmetric = false, bool fitTail = false, GFTail gft = GFTail()):
                symmetric(symmetric), beta(beta), fft(beta), deltaIt(static_cast<RealT>(MAX_T_BINS)/beta), fit(fitTail), tailFitted(false), gft(gft)
        {
            //tailFunc = defaultTail;//[](int n, int i) { return std::pow(1.0/ComplexT(0.,mFreqS(n, 10.)),i); 
            g_it = ImTG::Zero(MAX_T_BINS, SPINS);
            g_wn = MatG::Zero(MAX_M_FREQ, SPINS);
            g_it_std = ImTG::Zero(MAX_T_BINS, SPINS);
            g_wn_std = MatG::Zero(MAX_M_FREQ, SPINS);
            expansionCoeffs = {};
        }

        GreensFct(GreensFct const&) = delete; // TODO: implement

        virtual ~GreensFct(){}

        /*! Call this after setting \f$ G(\tau) \f$ 
         *  @return 1 if previously set, 0 otherwise
         */
        int markTSet(void){ fitTail(); return tSet++;}

        /*! Call this after setting \f$ G(i \omega_n) \f$.
         *  @return returns 1 if previously set, 0 otherwise
         */
        int markMSet(void){ fitTail(); return mSet++;}

        /*! Sets G(i \omega_n) = val
         *
         *  @param [in] n Matsubara frequency 
         *  @param [in] spin \f$ \sigma \in N \f$
         *  @param [in] val value of \f$ G_\sigma (\tau ) \f$
         */
        inline void setByMFreq(const int n, const int spin, const ComplexT val)
        {
            if(std::abs(n) >= g_wn.rows())
            {
                LOG(WARNING) << "out of bound for internal storage of Matsubara GF. Use tail instead!";
                return;
            }
            if(symmetric)
            {
                if(n<0) g_wn(-n,spin) = -val;
                else g_wn(n,spin) = val;
            }
            else
            {
                g_wn(n+g_wn.rows()/2,spin) = val;
            }
            tSet = 0;
            tailFitted = false;
        }


        /*! Sets G(i \omega_n) = val
         *
         *  @param [in] n Matsubara frequency 
         *  @param [in] spin \f$ \sigma \in N \f$
         *  @param [in] val value of \f$ G_\sigma (\tau ) \f$
         *  @param [in] std variance of data
         */
        inline void setByMFreq(const int n, const int spin, const ComplexT val, const ComplexT std)
        {
            if(symmetric)
            {
                if(n<0){
                    g_wn(-n,spin) = -val;
                    g_wn_std(-n,spin) = std;
                }
                else{
                    g_wn(n,spin) = val;
                    g_wn_std(n,spin) = std;
                }
            }
            else
            {
                g_wn(n+g_wn.rows()/2,spin) = val;
                g_wn_std(n+g_wn.rows()/2,spin) = std;
            }
            tSet = 0;
            tailFitted = false;
        }

        /*! Sets G(i \omega_n) = val
         *
         *  @param [in] g Matsubara Green's function (as MatG) 
         */
        void setByMFreq(const MatG &g)	{ g_wn = g; mSet = 1; tSet = 0; tailFitted = false;}

        /*! Sets G(i \omega_n) = val
         *
         *  @param [in] g Matsubara Green's function (as Eigen::ArrayXXd) 
         *  @param [in] std standard deviation of data 
         */
        void setByMFreq(const MatG &g, MatG &std)	{ g_wn = g; g_wn_std = std; mSet = 1; tSet = 0; tailFitted = false;}

        /*! Sets \f$ G(t) = val, G(0) = val \Rightarrow G(0^+) = val \f$
         *
         *  @param [in] t imaginary time \f$ \tau \f$
         *  @param [in] spin \f $\sigma \in N \f$
         *  @param [in] val value of Green's function for \sigma at \tau
         *  @param [in] std standard deviation of data
         */
        inline void setByT(const RealT t, const int spin, const RealT val, const RealT std = 0.0)
        {
            const int i = tIndex(t);
            //if(t==1) g_it(0,spin) = val;	// NOTE: redudant since tIndex shifts up from 0, remove?
            g_it(i,spin) = val;
            g_it_std(i,spin) = std;
            mSet = 0;
        }

        /*! Set G(\tau) = val
         *
         *  @param [in] g_it imaginary time Green's function (as ImTG)
         */
        void setByT(const ImTG &g) {g_it = g; markTSet(); tSet = 1; mSet = 0;}

        /*! Set G(\tau) = val
         *
         *  @param [in] g_it imaginary time Green's function (as ImTG)
         *  @param [in] std standard deviation of data
         */
        void setByT(const ImTG &g, ImTG &std) {g_it = g; g_it_std = std; markTSet(); tSet = 1; mSet = 0;}

        /*! get Matsubara Green's function \f$ G_\sigma(i \omega_n) \f$ by n.
         *  @param  n Matsubara frequency \f$ \omega_n \f$
         *  @param  spin \f$ \sigma \in N \f$
         *
         *  @return
         */
        inline ComplexT getByMFreq(const int n, const unsigned int spin) const
        {
            if(std::abs(n) >= g_wn.rows())
            {
                VLOG(5) << "using tail data instead of measured data, for wn: " << getMFreq(n) << " = " <<  getTail(getMFreq(n), spin);
                if(!tailFitted) LOG(WARNING) << "requesting tail without having called fit";
                return getTail(getMFreq(n), spin);
            }
            if(symmetric)
            {
                if(n<0) return -g_wn(-n,spin);
                else return g_wn(n,spin);
            }
            else
            {
                return g_wn(n+g_wn.rows()/2,spin);
            }
        }
        

        /*! Returns Matsubara frequency as function of the index n.
         */
        inline RealT getMFreq(const int n) const
        {
            return mFreqS(n, beta);
        }

        /*! @brief	imaginary time Green's function \f$ G_\sigma(\tau) \f$
         *  		storage order is t first (inner loop over t)
         *  
         *  @param [in] t imaginary time \f$ \tau \f$
         *  @param [in] spin \f$ \sigma \in N \f$
         *
         *  @return G_sigma(  \f$ \tau \f$ )
         */
        inline RealT operator() (const RealT t, const int spin) const
        {
            //REMARK: this implies that G(0) == G(0^-) 
            if(t == 0.0) return -1.0*g_it(g_it.rows()-1, spin);
            const int sgn = 2*(t>0)-1;
            //if(t<0){ t+=beta; return -1.0*g_it(tIndex(t), spin);}
            return sgn*g_it( tIndex( t + (t<0)*beta ), spin);
        }

        inline VectorT operator() (const VectorT t, const int spin) const
        {
            VectorT res(t.size());
            for(int i =0; i < t.size(); i++)
            {
                const int sgn = 2*(t[i]>0)-1;
                if(t[i] == 0.0) res(i) = -1.0*g_it(g_it.rows()-1, spin);
                else res[i] = sgn*g_it( tIndex( t[i] + (t[i]<0)*beta ), spin);
            }
            return res;
        }

        /*! get imaginary time Green's function using inperpolation.
         *
         *  @param [in] t imaginary time \f$ \tau \f$
         *  @param [in] spin \f$ \sigma \in N \f$
         *  @param [in] order of polynomial interpolation. TODO: inmpl. order > 1
         *
         *  @return G_sigma(  \f$ \tau \f$ )
         */
        RealT getByT(const RealT t, int spin, const int order) const
        {
            const int sgn = 2*(t>0)-1;
            const RealT tp = t + (t<0)*beta; 
            const int index0 = tIndex( tp );
            const RealT h = tp*static_cast<RealT>(MAX_T_BINS)/beta - index0;
            const int index1 = index0 < (MAX_T_BINS-1) ? (index0+1) : 1;
            const RealT val0 = sgn*g_it(index0, spin);
            const RealT val1 = sgn*g_it(index1, spin);
            //TODO: test!
            const RealT res = val0 + h*(val1 - val0);
            return res;
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
            res << std::fixed << std::setw(8)<< "\t mFreq \tSpin UP Re \t Spin UP Im \t Spin DOWN Re \t Spin DOWN Im \tSpin UP Re Err \t Spin UP Im Err \t Spin DOWN Re Err \t Spin DOWN Im Err \n";
            for(int n=0; n < MAX_M_FREQ; n++)
            {
                RealT wn = symmetric ? mFreqS(n, beta) : mFreq(n, beta);
                res << std::setprecision(8) << wn << "\t" <<  g_wn(n,UP).real() << "\t" << g_wn(n,UP).imag() << "\t" << g_wn(n,DOWN).real() << "\t" << g_wn(n,DOWN).imag() << "\t" <<  g_wn_std(n,UP).real() << "\t" << g_wn_std(n,UP).imag() << "\t" << g_wn_std(n,DOWN).real() << "\t" << g_wn_std(n,DOWN).imag() << "\n" ;
            }
            for(int n=MAX_M_FREQ; n < 2*MAX_M_FREQ*fit; n++)
            {
                RealT wn = symmetric ? mFreqS(n, beta) : mFreq(n, beta);
                ComplexT mExp[2];
                mExp[UP] = getTail( wn, UP);
                mExp[DOWN] = getTail( wn, DOWN);
                res << std::setprecision(6) << wn << "\t" <<  mExp[UP].real() << "\t" << mExp[UP].imag() << "\t" << mExp[DOWN].real() << "\t" << mExp[DOWN].imag() << "\n" ;
            }
            return res.str();
        }

        std::string getMaxEntString(int spin = UP) const
        {
            //if(!mSet) return "";
            std::stringstream res;
            res << std::fixed << std::setw(8);
            for(int n=0; n < MAX_M_FREQ; n++)
            {
                ComplexT err = (g_wn_std(n,spin) == ComplexT(0.0,0.0)) ? ComplexT(1.0/(2.0*(n*n+10.0)), 1.0/(2.0*(n*n+10.0))) : g_wn_std(n,spin);
                RealT wn = symmetric ? mFreqS(n, beta) : mFreq(n, beta);
                res << std::setprecision(8) << wn << "\t" << g_wn(n,spin).imag() << "\t" << err.imag() << "\n";
            }
            for(int n=MAX_M_FREQ; n < 5*MAX_M_FREQ*fit; n++)
            {
                RealT wn = symmetric ? mFreqS(n, beta) : mFreq(n, beta);
                ComplexT mExp[2];
                mExp[UP] = getTail( wn, spin);
                res << std::setprecision(8) << wn << "\t" << mExp[UP].imag() << "\t"<< 1.0/(MAX_M_FREQ*MAX_M_FREQ)<< "\n" ;
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
            //if(!tSet) return "";
            std::stringstream res;
            res << std::fixed << std::setw(8)<< "iTime \tSpin up \tSpin down \tSpin up Err\tSpin down Err\n";
            for(int i =0; i < MAX_T_BINS; i++)
            {
                res << std::setprecision(8) << std::setw(8) << (beta*i)/MAX_T_BINS << "\t"\
                    << g_it(i,UP)<< "\t" << g_it(i,DOWN) << "\t" << g_it_std(i,UP) << "\t" << g_it_std(i,DOWN) <<"\n";
            }
            return res.str();
        }

        /*! Does a FFT from Matsubara to imaginary time Green's function and stores it in internal storage.
        */
        void transformMtoT(void)
        {
            //if(symmetric) 
            transformMtoT(g_it);
            //if(fit)
            //   fitTail(); 
            //fft.transformMtoT_naive(this);
        }

        /*! Does a FFT from Matsubara to imaginary time Green's function and stores it in target.
         *  @param	[out]	Eigen::ArrayXXd where imaginary time Green's function will be stored
         */
        void transformMtoT(ImTG& target)
        {
            if(fit)
            {
                std::vector<std::array<RealT,2> > expC;
                fitTail();
                for(int i = 0; i < expansionCoeffs[0].size(); i++)
                {
                    expC.push_back({expansionCoeffs[0][i], expansionCoeffs[1][i]});
                }
                fft.transformMtoT(g_wn, target, expC, symmetric);
            }
            else
            {
                fft.transformMtoT(g_wn, target, symmetric);
            }
            markTSet();
        }

        /*! Does a FFT from imaginary time to Matsubara Green's function and stores it.
        */
        void transformTtoM(void){ fft.transformTtoM(g_it, g_wn); markMSet(); }   

        /*! Does a FFT from imaginary time to Matsubara Green's function and stores it in target.
         *  @param	[out]	Eigen::ArrayXXcd where Matsubara Green's function will be stored
         */
        void transformTtoM(MatG& target){ fft.transformTtoM(g_it, target); markMSet(); }

        /*! Compute \f[ ()G(i \omega_n)^{-1} - s)^{-1}  \f] and also transform to \f[ G(\tau) \f].
         *  @param  [in]  s shift
         */
        inline void shift(const RealT s) { g_wn = ((g_wn.cwiseInverse() - s).cwiseInverse()).eval() ; markMSet(); transformMtoT();}

        inline bool isSymmetric(void) const {return symmetric;}
        inline bool hasTail(void) {return fit;}

        /*! Compute this^-1 (i wn) - other^-1 (i wn), with tails if available.
         *
         *  @param [in] other MatG
         */
        MatG invAdd(GreensFct& other)
        {
            MatG res(1,1);
            fitTail();
            other.fitTail();
            if(fit && other.hasTail())
            {
                MatG res(MAX_T_BINS,2); 
                if(symmetric && other.isSymmetric())
                {
                    res.topRows(MAX_M_FREQ) = g_wn.cwiseInverse() - other.getMGF().cwiseInverse();
                    for(int n=MAX_M_FREQ;n<MAX_T_BINS;n++)
                    {
                        const RealT wn = mFreq(n,beta);
                        res(n,DOWN) = 1.0/getTail(wn, DOWN) - 1.0/other.getTail(wn, DOWN);
                        res(n,UP) = 1.0/getTail(wn, UP) - 1.0/other.getTail(wn, UP);
                    }
                }
                else
                {
                    int offset = static_cast<int>( (MAX_T_BINS - MAX_M_FREQ)/2.0 );
                    for(int n=0;n<offset;n++)
                    {
                        const RealT wn = mFreq(n,beta);
                        res(n,DOWN) = 1.0/getTail(wn, DOWN) - 1.0/other.getTail(wn, DOWN);
                        res(n,UP) = 1.0/getTail(wn, UP) - 1.0/other.getTail(wn, UP);
                    }
                    res.block(offset,0,MAX_M_FREQ,2) = g_wn.cwiseInverse() - other.getMGF().cwiseInverse();
                    for(int n=MAX_M_FREQ+offset;n<MAX_T_BINS;n++)
                    {
                        const RealT wn = mFreq(n,beta);
                        res(n,DOWN) = 1.0/getTail(wn, DOWN) - 1.0/other.getTail(wn, DOWN);
                        res(n,UP) = 1.0/getTail(wn, UP) - 1.0/other.getTail(wn, UP);
                    }
                }
            }
            else
            {
                MatG res = g_wn.cwiseInverse() - other.getMGF().cwiseInverse();
            }
            return res;
        }

        const MatG& getMGF() const { return g_wn; }
        const ImTG& getItGF() const { return g_it; }

        void setParaMagnetic(void)
        {
            g_wn.rightCols(1) = 0.5*(g_wn.leftCols(1) + g_wn.rightCols(1)).eval();
            g_wn.leftCols(1) = g_wn.rightCols(1).eval();
            fitTail();
            transformMtoT();
        }

        /*! Fits the expansion coefficients for the matsubara Green's function.
         *  This is used to interally improve accuracy.
         *
         *  @return void
         */
        void fitTail(void)
        {
            expansionCoeffs = {};
            //TODO: tailFitted flag not working correctly
            if(!fit) return;
            int nSamples = gft.last-gft.first + 1;
            Eigen::Matrix<ComplexT, Eigen::Dynamic, Eigen::Dynamic> A;
            Eigen::Matrix<ComplexT, Eigen::Dynamic, 1> rhs;
            Eigen::Matrix<RealT, 1, Eigen::Dynamic> fittedC;
            A.resize(nSamples, gft.nC);
            rhs.resize(nSamples);
            fittedC.resize(gft.nC);
            for(int f = 0; f < _CONFIG_spins; f++)
            {
                for(int i = 0; i < nSamples && (i + gft.first < _CONFIG_maxMatsFreq); i++)
                {
                    if(i + gft.first > _CONFIG_maxMatsFreq) LOG(ERROR) << "too many values for tail fit requested";
                    for(int j = 0; j < gft.nC; j++)
                    {
                        A(i,j) = gft.fitFct( gft.first + i , j, beta);
                        rhs(i) = getByMFreq( gft.first + i, f);
                    }
                }
                fittedC = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs).real().transpose();
                for(int j = 0; j < gft.nC; j++)
                {
                    if(!(isSymmetric() && j%2 == 0))
                        expansionCoeffs[f][j]  = fittedC(j);
                    else
                        expansionCoeffs[f][j]  = 0.;
                }
                VLOG(4) << "Tail coef f = " << f << ": " << expansionCoeffs[f][0] << ", " << expansionCoeffs[f][1] << ", " << expansionCoeffs[f][2] << ", " << expansionCoeffs[f][3] << ", " << expansionCoeffs[f][4];
            }
            tailFitted = true;
        }

        ComplexT getTail(const RealT wn, int spin) const
        {
            ComplexT iwn(0.0, wn);
            ComplexT res(0.,0.);
            for(int i = 0; i < expansionCoeffs[spin].size(); i++)
                res += expansionCoeffs[spin][i]/std::pow(iwn,i);
            return res;
        }

        RealT getTailCoef(const unsigned int i, int spin)
        {
            return expansionCoeffs[spin][i];
        }

        RealT getFtTail(const RealT tau, const int spin) const
        {
            RealT res = expansionCoeffs[spin][0];
            if(TAIL_ORDER > 1)
                res += -expansionCoeffs[spin][1]/2.0;
            if(TAIL_ORDER > 2)
                res +=  expansionCoeffs[spin][2]*(2.0*tau - beta)/4.0;
            if(TAIL_ORDER > 3)
                res += expansionCoeffs[spin][3]*tau*(beta-tau)/4.0;
            if(TAIL_ORDER > 4)
                res += expansionCoeffs[spin][4]*(2.0*tau - beta)*(2.0*tau*tau - 2.0*tau*beta - beta*beta)/48.0;
            return res;
        }

        protected:
        RealT beta;
        const bool symmetric;
        const bool fit;

        bool tailFitted;

        ImTG	g_it; // col major -> spin outer loop
        MatG	g_wn; // col major -> spin outer loop
        ImTG	g_it_std; // col major -> spin outer loop
        MatG	g_wn_std; // col major -> spin outer loop
        FFT 	fft;
        RealT   deltaIt;
        GFTail  gft;

        std::array<std::array<RealT,TAIL_ORDER>, _CONFIG_spins> expansionCoeffs;

        //Config& conf;

        int mSet = 0;	// track whether matsuabra or iTime GF has been set
        int tSet = 0;	// TODO implement

        /*! computes index from imaginary time. G(0) => G(0^+), G(beta) => G(beta^-)
         *  @param [in] t imaginary time \tau
         *
         *  @return index for lookup in g_it
         */
        inline int tIndex(const RealT t) const {
            const int res = static_cast<int>(t*deltaIt);
            return res; // TODO: test for +(res == 0);
        }


    };


    //using GreensFct = GreensFct;	// useless for now, may be useful for later changes to GreensFct
    //using GreensFct = GreensFct<_CONFIG_maxTBins, _CONFIG_maxMatsFreq, _CONFIG_spins>;
}	//end namespace DMFT

#endif
