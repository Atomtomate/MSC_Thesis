#ifndef INTEGRALS_HPP_
#define INTEGRALS_HPP_

#include "Config.hpp"
#include <gsl/gsl_integration.h>


namespace DMFT
{
// TODO: interfaces for DOS, e, integral
namespace examples
{
namespace integrals
{
	struct Args
    {
        Args(const RealT a, const RealT  e0, const RealT t, const RealT mu, const RealT U):
            a(a), e0(e0), t(t), mu(mu), U(U) {}
        const RealT a, e0, t, mu, U;
        RealT mf;
        ComplexT sImp;
    };
	inline double betheNer(double x, Args& args)
    {
        return std::real(std::sqrt(4.0*args.t*args.t-x*x)/\
            ((2.0*boost::math::constants::pi<double>()*args.t*args.t)*(std::complex<double>(args.mu - x,args.mf) - args.sImp)));
    }
	inline double betheNei(double x, Args& args)
    {
        return std::imag(std::sqrt(4.0*args.t*args.t-x*x)/\
            ((2.0*boost::math::constants::pi<double>()*args.t*args.t)*(std::complex<double>(args.mu - x,args.mf) - args.sImp)));
    }
    // t == p[0], mu == p[1], U == p[2], mf == p[3], sImp_re == p[4], sImp_imag == p[5]
	inline double GSL_betheNer(double x, void* args)
    {
        double* p = static_cast<double*>(args);
        return std::real(std::sqrt(4.0*p[0]*p[0]-x*x)/\
            (2.0*boost::math::constants::pi<double>()*p[0]*p[0]*std::complex<double>(p[1] - x - p[4],p[3] - p[5])));
    }
	inline double GSL_betheNei(double x, void* args)
    {
        double* p = static_cast<double*>(args);
        return std::imag(std::sqrt(4.0*p[0]*p[0]-x*x)/\
            (2.0*boost::math::constants::pi<double>()*p[0]*p[0]*std::complex<double>(p[1]  - x - p[4],p[3] - p[5])));
    }


    // Simple Cubic Lattice
	inline double sc2Dr(double x, double y, Args& args)
    {
        return std::real(1.0/(std::complex<double>(args.mu - args.e0  + \
            args.t*(std::cos(x*args.a) + std::cos(y*args.a)),args.mf) - args.sImp));
    }
	inline double sc2Di(double x, double y, Args& args)
    {
        return std::imag(1.0/(std::complex<double>(args.mu - args.e0  + \
            args.t*(std::cos(x*args.a) + std::cos(y*args.a)),args.mf) - args.sImp));
    }


    // t == p[0], mu == p[1], mf == p[2], sImp_re == p[3], sImp_imag == p[4], a == p[5] , e0 == p[6], xVal == p[7]
    inline double GSL_sc2Dr_h(double y, void* args){
        double* p = static_cast<double*>(args);
        return std::real(1.0/(std::complex<double>(p[1] + p[6]  - \
            p[0]*(std::cos(p[7]*p[5]) + std::cos(y*p[5])) ,p[2])));
    }
    inline double GSL_sc2Di_h(double y, void* args){
        double* p = static_cast<double*>(args);
        return std::imag(1.0/(std::complex<double>(p[1] + p[6]  - \
            p[0]*(std::cos(p[7]*p[5]) + std::cos(y*p[5])) ,p[2])));
    }
    // t == p[0], mu == p[1], mf == p[2], sImp_re == p[3], sImp_imag == p[4], a == p[5] , e0 == p[6], lim == p[7]
	inline double GSL_sc2Dr(double x, void* args)
    {

        gsl_function f;
        gsl_integration_cquad_workspace *ws;
        if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) {
            LOG(ERROR) << "main: call to gsl_integration_cquad_workspace_alloc failed.";
        }
        double res;
        f.function = GSL_sc2Dr_h;
        double* p = static_cast<double*>(args);
        double sub_p[8] = {p[0],p[1],p[2],p[3],p[4],p[5],p[6],x};
        f.params = sub_p;
        if ( gsl_integration_cquad( &f,-p[7],p[7],1.0e-10,1.0e-10,ws,&res,NULL,NULL) != 0 ) {
            LOG(ERROR) <<  "main: call to gsl_integration_cquad in GSL_sc2Dr failed.";
        }
        gsl_integration_cquad_workspace_free( ws );
        return res;
    }
	inline double GSL_sc2Di(double x, void* args)
    {
        gsl_function f;
        gsl_integration_cquad_workspace *ws;
        if ( ( ws = gsl_integration_cquad_workspace_alloc( 200 ) ) == NULL ) {
            LOG(ERROR) << "main: call to gsl_integration_cquad_workspace_alloc failed.";
        }
        double res;
        f.function = GSL_sc2Di_h;
        double* p = static_cast<double*>(args);
        double sub_p[8] = {p[0],p[1],p[2],p[3],p[4],p[5],p[6],x};
        f.params = sub_p;
        if ( gsl_integration_cquad( &f,-p[7],p[7],1.0e-10,1.0e-10,ws,&res,NULL,NULL) != 0 ) {
            LOG(ERROR) <<  "main: call to gsl_integration_qagi in GSL_sc2Di failed.";
        }
        gsl_integration_cquad_workspace_free( ws );
        return res;
    }



	inline double sc3Dr(double x, double y, double z, Args& args)
    {
        return std::real(1.0/(std::complex<double>(args.mu - args.e0  + \
            2.0*args.t*(std::cos(x*args.a) + std::cos(y*args.a) + std::cos(z*args.a)),args.mf) - args.sImp));
    }
	inline double sc3Di(double x, double y, double z, Args& args)
    {
        return std::imag(1.0/(std::complex<double>(args.mu - args.e0  + \
            2.0*args.t*(std::cos(x*args.a) + std::cos(y*args.a) + std::cos(z*args.a)),args.mf) - args.sImp));
    }
}   // end namespace integrals
	namespace DOS{
		inline double Bethe_Elliptical(const double t, const double e){
			const double t2 = t*t;
			return std::sqrt(4.0*t2-e*e)/(2.0*boost::math::constants::pi<double>()*t2);
		}
        inline std::complex<double> Bethe_Elliptical_Hilbert(const double t2, const std::complex<double> zeta)
        {
            const int s = 2*(std::imag(zeta) > 0) - 1;
            return (zeta - static_cast<double>(s)*std::sqrt(zeta*zeta - 4.0*t2))/(2.0*t2);
        }

        inline double SC(const double t, const double e){
            return 1.0/(t*std::sqrt(2.0*boost::math::constants::pi<double>()))*std::exp(-e*e/(2.0*t));
        }
	}

	namespace energies{
		template<int D>
		inline double scND(const std::array<double,D> k,const double e0, const double t, const double a){
			double result = 0.0;
			for(int i=0;i<D;i++) result += std::cos(k[i]*a);
			return e0 + 2.0*t*result;
		}
	}
} //end namespace examples

} //end namespace DMFT


#endif
