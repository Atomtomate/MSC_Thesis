#include "Memoization.hpp"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <boost/timer/timer.hpp> // timers
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"

typedef std::complex<double> cdT;
using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;
/*! compile with icpc/clang/g++ -std=c++14 test_Memoization.cpp -I/home/stobbe/software/include/boost -L/home/stobbe/software/lib -lboost_timer -lboost_system
 */

inline cdT iexp(double a, double phi){
	return std::exp(a + std::complex<double>(0.0,phi));
}

void test_perf(){
	// setting up primitive container as comparison
	int N = 2000;
	cdT **test_arr = new cdT*[N];
	for(int i = 0; i < N; ++i) {
    	test_arr[i] = new cdT[N];
	}
	boost::iostreams::stream< boost::iostreams::null_sink > nullStream( ( boost::iostreams::null_sink() ) );

	// setting up timer
	nanosecond_type last(0);
	cpu_timer timer;

	// setting up memoize
	auto memo_iexp = Memoize::memoize( iexp );

	for(int a = 0; a < N; a++){
		for(int phi_i = 0; phi_i < N; phi_i++ ){
			double phi = phi_i*2*3.14159265/N;
			test_arr[a][phi_i] = iexp(a,phi);
				
		}
	}

  	cpu_times elapsed_times(timer.elapsed());
  	nanosecond_type elapsed(elapsed_times.system + elapsed_times.user);	
	std::cout << "naive array initialized, time:\t" << timer.format(5) << std::endl;
	timer.start();
	for(int a = 0; a < N; a++){
		for(int phi_i = 0; phi_i < N; phi_i++ ){
			double phi = phi_i*2*3.14159265/N;
			memo_iexp(a,phi);
				
		}
	}
	std::cout << "memoized initialized, time:\t" << timer.format(5) << std::endl;
	std::srand(std::time(0));
	timer.start();

	for(int phi_i = 0; phi_i < N; phi_i++ ){
		for(int a = 0; a < N; a++){
			auto a_r = std::rand()%N;
			auto phi_r = std::rand()%N;
			double phi = phi_r*2*3.14159265/N;
			nullStream << test_arr[a_r][phi_r] << phi;
				
		}
	}
	std::cout << "naive array read, time:\t" << timer.format(5) << std::endl;
	timer.start();

	for(int phi_i = 0; phi_i < N; phi_i++ ){
		for(int a = 0; a < N; a++){
			auto a_r = std::rand()%N;
			auto phi_r = std::rand()%N;
			double phi = phi_r*2*3.14159265/N;
			nullStream << memo_iexp(a_r,phi) << phi;
				
		}
	}
	std::cout << "memoized read, time:\t" << timer.format(5) << std::endl;


	for(int i = 0; i < N; ++i) {
    	delete [] test_arr[i];
	}
	delete [] test_arr;	
}
//iexp_hash()

void test_umap_perf(){
	// setting up primitive container as comparison
	int N = 2000;
	cdT **test_arr = new cdT*[N];
	for(int i = 0; i < N; ++i) {
    	test_arr[i] = new cdT[N];
	}
	boost::iostreams::stream< boost::iostreams::null_sink > nullStream( ( boost::iostreams::null_sink() ) );

	// setting up timer
	nanosecond_type last(0);
	cpu_timer timer;

	// setting up memoize
	auto memo_iexp = Memoize::memoize( iexp );

	for(int a = 0; a < N; a++){
		for(int phi_i = 0; phi_i < N; phi_i++ ){
			double phi = phi_i*2*3.14159265/N;
			test_arr[a][phi_i] = iexp(a,phi);
				
		}
	}

  	cpu_times elapsed_times(timer.elapsed());
  	nanosecond_type elapsed(elapsed_times.system + elapsed_times.user);	
	std::cout << "naive array initialized, time:\t" << timer.format(5) << std::endl;
	timer.start();
	for(int a = 0; a < N; a++){
		for(int phi_i = 0; phi_i < N; phi_i++ ){
			double phi = phi_i*2*3.14159265/N;
			memo_iexp(a,phi);
				
		}
	}
	std::cout << "memoized initialized, time:\t" << timer.format(5) << std::endl;
	std::srand(std::time(0));
	timer.start();

	for(int phi_i = 0; phi_i < N; phi_i++ ){
		for(int a = 0; a < N; a++){
			auto a_r = std::rand()%N;
			auto phi_r = std::rand()%N;
			double phi = phi_r*2*3.14159265/N;
			nullStream << test_arr[a_r][phi_r] << phi;
				
		}
	}
	std::cout << "naive array read, time:\t" << timer.format(5) << std::endl;
	timer.start();

	for(int phi_i = 0; phi_i < N; phi_i++ ){
		for(int a = 0; a < N; a++){
			auto a_r = std::rand()%N;
			auto phi_r = std::rand()%N;
			double phi = phi_r*2*3.14159265/N;
			nullStream << memo_iexp(a_r,phi) << phi;
				
		}
	}
	std::cout << "memoized read, time:\t" << timer.format(5) << std::endl;


	for(int i = 0; i < N; ++i) {
    	delete [] test_arr[i];
	}
	delete [] test_arr;	
}

int main(void)
{	
	//test_perf();		
	test_umap_perf();
    return 0;
}
