#ifndef DMFT_SOLVER_HPP_
#define DMFT_SOLVER_HPP_
#include "GreensFct.hpp"
#include "WeakCoupling.hpp"
#include "KInt.hpp"
#include "IOhelper.hpp"
#include "ImpSolver.hpp"
#include "integrals.hpp"

#include <string>
#include <gsl/gsl_integration.h>

#define USE_GSL 1
//#define BETHE_LATTICE 1

namespace DMFT
{
    struct Params
    {
        Params(const RealT U,const RealT beta, RealT mu,const RealT e0):
            beta(beta), U(U), mu(mu), e0(e0) {}
        const RealT beta;		// inverse temperature
        const RealT	U;			// interaction strength
        RealT	mu;			// chemical potential
        const RealT	e0;			// ground state (for now half filling: U/2) 
    };
    /*!
     *  - guess \f$ G_0(\tau) \f$ 
     *  do {
     *		- compute \f$ \Sigma_{Imp}(i \omega_n) \f$ from \f$ G_0(\tau) \f$ 
     *	   		-# sample from Solver.update()
     *	   		-# extract \f$ G_{Imp}(i \omega_n) \f$ 
     *	   		-# \f$ \Sigma_{Imp}(i \omega_n) = G_0(i \omega_n)^{-1} - G_{Imp}(i \omega_n)^{-1} \f$ 
     *		- compute \f$ G_{Loc} (i \omega_n) \f$ 
     *		    - \f$ \int\limits_{\pi/a}^{\pi/a} [i \omega_n + \mu - \epsilon_k - \Sigma_{Imp} (i \omega_n)]^{-1} dk \f$ 
     *		- \f$ G_0(\tau) \f$ 
     *	    	- \f$ G_0(i \omega_n)^{-1} = G_{Loc} (i \omega_n)^{-1} + \Sigma_{Imp}(i \omega_n) \f$ 
     *	    	- FFT \f$ G(i \omega_n) \f$ to \f$ G(\tau)\f$ 
     *
     *  } while( G_loc not converged )
     *
     *	- issues:
     *		- to avoid calls to virtual functions we don't use polymorphism for now
     */
    template<class ImpSolver>
        class DMFTSolver
        {
            public:
                DMFTSolver(const std::string& outDir, const Config& config, RealT mixing, RealT D, ImpSolver &solver, GreensFct &G0, GreensFct &GImp, GreensFct &GLoc, const Params& p):
                    par(p), mixing(mixing),t(t), iSolver(solver),\
                    g0(G0),g0Info("G0"),gImp(GImp),gImpInfo("GImp"),gLoc(GLoc),gLocInfo("GLoc"),sImpGF(config.beta),sImpInfo("sImp"),gDelta(config.beta),gDeltaInfo("Delta"),\
                    ioh(outDir), kInt()
            {
                std::string tmp("G0_Guess");
                ioh.writeToFile(g0,tmp);
                ioh.addGF(g0,g0Info);
                ioh.addGF(gImp,gImpInfo);
                ioh.addGF(gLoc,gLocInfo);
                ioh.addGF(gDelta,gDeltaInfo);
                ioh.addGF(sImpGF,sImpInfo);
            };

                virtual ~DMFTSolver()
                {
#ifdef USE_GSL
                    gsl_integration_cquad_workspace_free( ws );
#endif
                }

                void writeToFile(void)
                {
                    LOG(INFO) << "Writing results";
                    ioh.writeToFile();
                }

                //TODO: use Math. Comp. 80 (2011), 1745-1767 for hilbert transform

                void update(const unsigned int iterations = 1, const unsigned long long updates = 500000)
                {
                    const double a_t1 = 1.0, e0 = 0.0;
                    examples::integrals::Args fArg(a_t1,e0,t,par.mu,par.U);

#ifdef USE_GSL
                    ws = NULL;
                    /* Initialize the workspace. */
                    if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) {
                        LOG(ERROR) << "main: call to gsl_integration_cquad_workspace_alloc failed.";
                    }
                    f_re.params = NULL;
                    f_im.params = NULL;
#ifdef BETHE_LATTICE
                    f_re.function = examples::integrals::GSL_betheNer;
                    f_im.function = examples::integrals::GSL_betheNei;
#else
                    f_re.function = examples::integrals::GSL_sc2Dr;
                    f_im.function = examples::integrals::GSL_sc2Di;
#endif
#endif


                    for(unsigned int dmftIt = 0; dmftIt < iterations; dmftIt++){
                        ioh.setIteration(dmftIt);

                        //TODO get G0 and gImp from solver
                        for(long unsigned int i=0; i <= 50; i++){
                            iSolver.update(updates/50.0);
                            LOG(INFO) << (2*i) << "%, rows: " << iSolver.expansionOrder();
                        }
                        LOG(INFO) << "finished sampling";
                        LOG(INFO) << "measuring impurity Greens function";
                        //iSolver.computeImpGF();
                        iSolver.computeImpGF();
                        LOG(INFO) << "Computing self energy from old Weiss GF and sampled impurity GF";
                        //REMARK: DMFTSolver is friend of GF to allow this to be inlined. encapsulate this?

                        MatG sImp = (g0.getMGF().cwiseInverse() - gImp.getMGF().cwiseInverse()) + par.U/2.0;
                        sImpGF.setByMFreq(sImp);
                        LOG(INFO) << "Computing new Weiss Green's function";
                        const RealT lim = boost::math::constants::pi<RealT>()/a_t1;       // SC
                        for(int s=0;s<_CONFIG_spins;s++){
                            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                                RealT mf = mFreq(n,par.beta);
#ifdef USE_GSL            
                                //TODO: move GSL integration to KInt
                                double res_re, abserr_re,res_im, abserr_im;
                                size_t neval_re, neval_im;
#ifdef BETHE_LATTICE
                                double gsl_args[6] = {t, par.mu, par.U, mf,std::real(sImp(n,s)), std::imag(sImp(n,s)) };
                                RealT lim = 2.0*t;       // bethe bandwidth
#else
                                double gsl_args[8] = {t, par.mu, mf,std::real(sImp(n,s)), std::imag(sImp(n,s)),1.0,0.0, lim};
#endif
                                f_re.params = gsl_args;
                                f_im.params = gsl_args;
                                if ( gsl_integration_cquad( &f_re, -lim,lim,1.0e-10,1.0e-10,ws,&res_re,NULL,NULL) != 0 ) {
                                    LOG(ERROR) <<  "main: call to gsl_integration_cquad with f_re failed.";
                                }
                                if ( gsl_integration_cquad( &f_im, -lim,lim,1.0e-10,1.0e-10,ws,&res_im,NULL,NULL) != 0 ) {
                                    LOG(ERROR) <<  "main: call to gsl_integration_cquad with f_im failed.";
                                }
                                gLoc.setByMFreq(n,s, ComplexT(res_re, res_im));
#else
#ifdef BETHE_LATTICE
                                fArg.mf = mf;
                                fArg.sImp= sImp(n,s);
                                auto res_R = kInt.integrateGL<1024, RealT, examples::integrals::Args&>(
                                        examples::integrals::betheNer,-2.0*fArg.t ,2.0*fArg.t, fArg );
                                auto res_I = kInt.integrateGL<1024, RealT, examples::integrals::Args&>(
                                        examples::integrals::betheNei,-2.0*fArg.t ,2.0*fArg.t ,fArg );
                                gLoc.setByMFreq(n,s, ComplexT(res_R, res_I));
#else
                                int kPoints = 500;
                                ComplexT result(0.0,0.0);
                                for(double kx = -lim, kx < lim, kx+= (2.0*lim)/kPoints)
                                {
                                    for(double ky = -lim, ky < lim, ky+= (2.0*lim)/kPoints)
                                    {
                                        result += 1.0/(ComplexT(par.mu - 2.0*(std::cos(kx) + std::cos(ky)),mf) );
                                    }
                                }
                                result = result/static_cast<RealT>(kPoints*kPoints);
                                gLoc.setByMFreq(n,s, result);
#endif
#endif
                                gDelta.setByMFreq(n,s, ComplexT(par.mu, mf) - 1.0/gLoc.getByMFreq(n,s) );
                                // TODO: BZ volume?
                            }
                        }

                        sImpGF.transformMtoT();
                        gLoc.markMSet();
                        gLoc.transformMtoT();
                        gDelta.markMSet();
                        gDelta.transformMtoT();

                        LOG(INFO) << "Writing results";
                        ioh.writeToFile();
                        //TODO: overload =Operator and use expression templates
                        //MatG tmp = mixing*g0.getMGF() + (1-mixing)*((gLoc.getMGF().cwiseInverse() + sImp).cwiseInverse());
                        MatG tmp = (gLoc.getMGF().cwiseInverse() + sImp).cwiseInverse();

                        for(int s=0;s<_CONFIG_spins;s++){
                            for(int n=0;n<_CONFIG_maxMatsFreq;n++){
                                RealT mf = mFreq(n,par.beta);
                                g0.setByMFreq(n,s, ComplexT(par.mu - par.U/2.0,mf) - (t/2.0)*(t/2.0)*gImp.getByMFreq(n,s) );
                            }
                        }
                        g0.markMSet();
                        //TODO: impose paramagnetic (gup+gdown) /2
                        LOG(INFO) << "Transforming G0 from Matsubara to imaginary time";

                        g0.transformMtoT();

                    }
                }


            private:
                const Params&	par;
                IOhelper	ioh;
                RealT mixing;
                RealT t;

                ImpSolver& iSolver;
                LogInfos	g0Info;
                LogInfos	gImpInfo;
                LogInfos	gLocInfo;
                LogInfos	sImpInfo;
                LogInfos    gDeltaInfo;
                GreensFct&	g0;
                GreensFct&	gImp;
                GreensFct&	gLoc;
                GreensFct  sImpGF;
                GreensFct   gDelta;
                utility::KInt		kInt;

#ifdef USE_GSL
                gsl_function f_re;
                gsl_function f_im;
                gsl_integration_cquad_workspace *ws;
#endif

        };
} //end namespace DMFT

//template classes need to be implemented inside header



/*ComplexT res = 0;
  const int kPoints = 100;
  for(RealT k_y=-PI/a_t1;k_y<PI/a_t1;k_y+= (2.0*PI/a_t1)/kPoints )
  for(RealT k_x=-PI/a_t1;k_x<PI/a_t1;k_x+= (2.0*PI/a_t1)/kPoints ){
//sc
res += 1.0/(ComplexT(par.mu - e0  + 2.0*t_t1*(std::cos(k_x*a_t1) + std::cos(k_y*a_t1)),mf) - sImp(n,s));
}*/
#endif
