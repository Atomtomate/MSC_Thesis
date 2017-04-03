#include "KInt.hpp"

#include <boost/array.hpp>
//#include <boost/numeric/odeint.hpp>
#include <boost/timer/timer.hpp>
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"
#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;

using namespace DMFT;

namespace utility
{

	namespace examples{
		namespace DOS{
			inline double elliptical(const double t, const double e){
				const double t2 = t*t;
				return std::sqrt(4.0*t2-e*e)/(2.0*3.14159265*t2);
			}
		}

		namespace energies{
			template<int D>
			double scND(const std::array<double,D> k,const double e0, const double t, const double a){
				double result = 0.0;
				for(int i=0;i<D;i++) result += std::cos(k[i]*a);
				return e0 + 2.0*t*result;
			}




		//	inline double bcc3d
		//   --:		return e0 + 8.0*t_t1*std::cos(k_x*a_t1/2.0)*std::cos(k_y*a_t1/2.0)*std::cos(k_z*a_t1/2.0);
		// fcc:  e0 + 4.0*t_t1*std::cos(k_x*a_t1/2.0)*std::cos(k_y*a_t1/2.0)+ std::cos(k_x*a_t1/2.0)*std::cos(k_z*a_t1/2.0)+ std::cos(k_y*a_t1/2.0)*std::cos(k_z*a_t1/2.0)
		}

		namespace integrals{
			void test_cos(const std::array<double,1> x, std::array<double,1>& dx, double t){
				dx[0] = std::cos(t);
			}

			/*void rhs(const std::array<double,1> x, std::array<double,1> &dx){
				dx[0] = 1.0/(energies::sc1d(e0,x[0],t,a) - c);
			}*/



		}
	}

// compile with:clang -std=c++14 test_KInt.cpp -I../includes -I/home/stobbe/software/include/boost -I../libs -I../libs/sprng5/include -L../libs/sprng5/lib -L/home/stobbe/software/lib -lboost_timer -lboost_system -lstdc++ -lm


//TODO:	mean over >1000 VERY large OR mean over >MAXINT number of samples numbers (check for overflow)
//TODO: test against GSL integration routines
inline RealT f1(double x, double&& a)
{
	return a+std::sin(x); //std::sin(x[0]);
}

inline RealT f2(double x, double y, double&& a)
{
	return 1.0/(a + 2*(std::cos(x) + std::cos(y)));
}
inline RealT f3(double x, double y, double z, double&& a)
{
	return 1.0/(2*a + 2*(std::cos(x) + std::cos(y) + std::cos(z)));
}
void printTest(std::array< double , 1 > &x, const double t){
	std::cout << t << "\t" << x[0] << "\t" << std::cos(t) << std::endl;
}



}// end namespace utilty

using namespace utility;

/*void test_2D_GL(void)
{

}*/

int main( int argc, char** argv ){
	// setting up environment
	double s = 0.0;
	double si = 4.0/10000.0;
	boost::iostreams::stream< boost::iostreams::null_sink > nullStream( ( boost::iostreams::null_sink() ) );

	double params[] = {-2.3,1.2,1.5,-13.4};		//e0,t,a,c
	std::array<RealT, 4> paramsArr = {-2.3,1.2,1.5,-13.4};


	// setting up timer
	nanosecond_type last(0);
	cpu_timer timer;


	s = 4.0*s/10000.0;
	std::cout << "1D direct sum: " << s << std::endl;
	std::cout << "using WeightedNDAcc now" << std::endl;
	// === 1D test ===
	// Riemann acc
	WeightedNDAcc<1> k1;
	std::array<RealT,1> min1 = {0};
	std::array<RealT,1> max1 = {4};
	std::array<unsigned long,1> N1 = {100};
  	cpu_times elapsed_times(timer.elapsed());
  	nanosecond_type elapsed(elapsed_times.system + elapsed_times.user);
	std::cout << "1D weighted riemann sum: " \
		 << ". Time: " <<  timer.format(5) << std::endl;
	timer.start();

	// GL acc
	KInt ki;
	double tmp = 4.0;
	std::cout << "Integrate 4+sin(x) from 0 to 10: " << ki.integrateGL<6>(f1,0.0,10.0,tmp ) << std::endl;
	std::cout << "Integrate 4+sin(x) from 0 to 10: " << ki.integrateGL<7>(f1,0.0,10.0,tmp ) << std::endl;

	tmp = 5.0;
	std::cout << "Integrate 1/(5 + 2*(cos(x) + cos(y))), x-[-10,10], y-[-6,6] (expected: 64.5456): ";
	std::cout << ki.integrateGL<128>(f2,-10.0,10.0,-6.0,6.0,tmp ) << std::endl;


	std::cout << "Integrate 1/(2*5 + 2*(cos(x) + cos(y) + cos(z))), x-[-10,10], y-[-1,1], z-[-2,5] (expected: 24.9884): " <<
		ki.integrateGL<32>(f3,-10.0,10.0,-1.0,1.0,-2.0,5.0,tmp) << std::endl;
	// === 2D test ===
/*	WeightedNDAcc<2> k2;
	std::array<RealT,2> min2 = {0,0};
	std::array<RealT,2> max2 = {4,9};
	std::array<unsigned long,2> N2 = {100,100};
	std::cout << "2D weighted riemann sum: " <<k2.sumKPoints(f2, min2, max2, N2) << ". Time: " <<  timer.format(5) <<std::endl;
	timer.start();

	double t0 = 0.0;
	std::array<double,1> x = {0.0};
	//boost::numeric::odeint::integrate(examples::integrals::test_cos, x, t0, 3.14159265/2.0, 0.0001, printTest);
	std::cout << x[0] << std::endl;*/
	return 0;
}
