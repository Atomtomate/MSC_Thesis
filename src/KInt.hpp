#ifndef KINT_H_
#define KINT_H_
#include "Config.hpp"
#include "GLWeights.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>
#include <boost/accumulators/framework/parameters/weight.hpp>
#include <boost/accumulators/framework/accumulators/external_accumulator.hpp>
#include <boost/accumulators/statistics/weighted_sum.hpp>
//#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <iostream>


/*!	This class provides methods for integration over momenta or energies
 *	
 *	Problem:
 *		Parameters are: 
 *		- a: lattice spacing
 *		- \f$ D(\epsilon) \f$: density of states (DOS)
 *		- Momentum, energy dependent function. In most cases this will be the impurity Green's \f$ G(k,i\omega_n) \f$
 *		- Some constant offset. In most cases \f$ i \omega_n + \mu - \Simga(i \omega_n)\f$
 *		- function and limits for the summand e.g.: 
 *		\f[
 *			\int\limits_{-\frac{\pi}{a_1}}^\frac{\pi}{a_1}  \cdots \int\limits_{-\frac{\pi}{a_D}}^\frac{\pi}{a_D}  \frac{1}{i \omega_n - \epsilon_k + \mu - \Sigma (i \omega_n)} d k_D \cdots dk_1
 *		\f]
 *		or
 *		\f[
 *			\int\limits_{-\infty}^\infty  \frac{D(\epsilon)}{i \omega_n - \epsilon + \mu - \Sigma (i \omega_n)} d \epsilon 
 *		\f]
 *
 *	Methods:
 *
 *	- summation (this constructs some general weighted D dimensional sum)
 *		- Riemann sum
 *		- Gauss-Legendre
 *		
 *	- not implemented yet:
 *	- reformulation to ODE and solution via boost ODE int
 * 
 *	- Linear tetrahedral method for 1,2,3 D
 * 		- issues: can break symmetry, gamma point must be included
 *		- use MFEM?
 *		.
 *
 *	- TODO: let user decide wether to use kahan summation or not
 *	- TODO: provide data structure for additional parameters
 */
namespace utility{
//using KahanTag = boost::accumulators::tag::weighted_sum_kahan;
//using KahanSumAccT = boost::accumulators::stats<KahanTag>;
using TagType = boost::accumulators::tag::weighted_sum;
using SumAccT = boost::accumulators::stats<TagType>;
template <typename T>
using AccT = boost::accumulators::accumulator_set<T, SumAccT, T >;


namespace detail
{
	/*!	this recursively constructs the depth of the nested for loops
	 *	for the innermost loop the partially specialized struct Internal<D,0> is called
	 */
	template<unsigned int D, unsigned int ND, typename T, typename RetT>
	struct Internal
	{
		static void sumKPoints(RetT (*summand)(std::array<T,D> x),const std::array<T,D> &min,const std::array<T,D> &incs,\
			 const std::array<unsigned long, D> &N,std::array<T, D> &xVec, AccT<RetT> &acc, const RetT &weight)
		{
			RetT xi=min[ND];
			for(unsigned int n=0; n<N[ND];n++)
			{
				xVec[ND] = xi;
				Internal<D,ND-1, T, RetT>::sumKPoints(summand, min, incs, N, xVec, acc, weight);
				xi+=incs[ND];
			}
		}
	};

	template<unsigned int D, typename T, typename RetT>
	struct Internal<D,0,T,RetT>
	{
		static void sumKPoints(RetT (*summand)(std::array<T,D> x),const std::array<T,D> &min,const std::array<T,D> &incs,\
			 const std::array<unsigned long, D> &N, std::array<T, D> &xVec, AccT<RetT> &acc, const RetT &weight)
		{
			RetT xi=min[0];
			for(unsigned int n=0; n<N[0];n++)
			{
				xVec[0] = xi;
				acc(summand(xVec), boost::accumulators::weight = weight);
				xi+=incs[0];
			}
		}
	};
}


//TODO: let user specify weights, TODO: let user choose between kahan and normal sum

template<unsigned int D>
class WeightedNDAcc
{
	public:
		/*!	@brief	instantiates an accumulator for D nested for loops \f$ \sum\limits_{min_1}^{max_1} ... \sum\limits_{min_D}^{max_D} f \f$
		 *			with \f$ f(T x_1, ... ,T x_D) -> T \f$
		 *
		 *	@param	summand		function pointer over which to sum
		 *	@param	min			array of start values (one element for each dimension)
		 *	@param	max			array of values for the upper limit (one element for each dimension)
		 *	@param	N			array with number of steps (one element for each dimension)
		 *
		 *	@return	accumulated value
		 */
		template<class RetT, class T >
		RetT sumKPoints(RetT (*summand)(std::array<T, D> x), std::array<T, D> &min, std::array<T, D> &incs, std::array<T, D> &N, std::array<T, D> &xVec) const
		{
			//acc<summand type, method, weight type>
			AccT<RetT> acc;
			RetT weight = 1.0;
			for(auto el : incs) weight *= el;

			detail::Internal<D,D-1,T,RetT>::sumKPoints(summand, min, incs, N, xVec, acc, weight);
			return boost::accumulators::sum(acc);
		}

};

/*! This class provides some integration Methods for up to 3 dimensions. 
 *	For more than 3 dimensions one should resort to MC methods (thay are provided by GSL)
 *	TODO: this can somehow be done with std::bind and function pointers...
 *	TODO: dynamic cache for GLWeights
 *  TODO: better matching for arguments
 *	as of now there is a lot of redundant code.
 */
class KInt
{
	public:
		KInt() {}
		/*template<unsigned order, typename RetT, typename... Args>
		RetT integrateGL(RetT (*f)(double , double, Args&&...), double from_x, double to_x, double from_y, double to_y, Args... args)
		{
			auto g = boost::bind(f,std::placeholders::_1,std::placeholders::_2, args...);
			intGL<order, RetT>(g,from_x,to_x,from_y,to_y);
		}

		template<unsigned order, typename RetT, typename... Args>
		RetT integrateGL(RetT (*f)(double x, Args&&... args), double from_x, double to_x, Args... args)
		{
			boost::function<RetT(double)> fb = std::bind(f,_1,args...);
			intGL<order, RetT>(fb,from_x,to_x);
		}*/



		template<unsigned order, typename RetT, typename... Args>
		RetT integrateGL(RetT (*f)(double x, double y, double z, Args&&... ),\
			double from_x, double to_x, double from_y, double to_y, double from_z, double to_z, Args... args)
		{
			AccT<RetT> acc;
			GLWeights::GLWeights<order> w;
			if(! GLWeights::defined<GLWeights::GLWeights<order>>::value)
			{
				//LOG(ERROR) << "no Legendre polynomial cache for this order!";
				std::cout << "error, order not cached" << std::endl;
			}
			const double c1x = (to_x-from_x)/2.0;
			const double c2x = (to_x+from_x)/2.0;
			const double c1y = (to_y-from_y)/2.0;
			const double c2y = (to_y+from_y)/2.0;
			const double c1z = (to_z-from_z)/2.0;
			const double c2z = (to_z+from_z)/2.0;
			int i = 0, j = 0, k = 0;
			const unsigned int limit = (order+1)/2;
			if(order%2)
			{
				std::cout << "error, uneven order not implemented" << std::endl;
			}
			while(i<limit)
			{
				j = order%2;
				while(j < limit)
				{
					k = order%2;
					while(k < limit)
					{
						acc(f( c1x*w.xi[i] + c2x, c1y*w.xi[j] + c2y, c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						acc(f( c1x*w.xi[i] + c2x, c1y*w.xi[j] + c2y,-c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						acc(f( c1x*w.xi[i] + c2x,-c1y*w.xi[j] + c2y, c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						acc(f( c1x*w.xi[i] + c2x,-c1y*w.xi[j] + c2y,-c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						acc(f(-c1x*w.xi[i] + c2x, c1y*w.xi[j] + c2y, c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						acc(f(-c1x*w.xi[i] + c2x, c1y*w.xi[j] + c2y,-c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						acc(f(-c1x*w.xi[i] + c2x,-c1y*w.xi[j] + c2y, c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						acc(f(-c1x*w.xi[i] + c2x,-c1y*w.xi[j] + c2y,-c1z*w.xi[k] + c2z, std::forward<Args>(args)...),\
							boost::accumulators::weight = w.wi[i]*w.wi[j]*w.wi[k]);
						k += 1;
					}
					j += 1;
				}
				i += 1;
			}
			return c1x*c1y*c1z*boost::accumulators::sum(acc);
		}

		template<unsigned order, typename RetT, typename... Args>
		RetT integrateGL(RetT (*f)(double x, double y, Args&&... ),\
			double from_x, double to_x, double from_y, double to_y, Args... args)
		{
			AccT<RetT> acc;
			GLWeights::GLWeights<order> w;
			if(! GLWeights::defined<GLWeights::GLWeights<order>>::value)
			{
				//LOG(ERROR) << "no Legendre polynomial cache for this order!";
				std::cout << "error, order not cached" << std::endl;
			}
			const double c1x = (to_x-from_x)/2.0;
			const double c2x = (to_x+from_x)/2.0;
			const double c1y = (to_y-from_y)/2.0;
			const double c2y = (to_y+from_y)/2.0;
			int i = 0, j = 0;
			const unsigned int limit = (order+1)/2;
			if(order%2)
			{
				std::cout << "error, uneven order not implemented" << std::endl;
			}
			while(i<limit)
			{
				j = order%2;
				while(j < limit)
				{
					acc(f( c1x*w.xi[i] + c2x, c1y*w.xi[j] + c2y, std::forward<Args>(args)...),\
						boost::accumulators::weight = w.wi[i]*w.wi[j]);
					acc(f( c1x*w.xi[i] + c2x,-c1y*w.xi[j] + c2y, std::forward<Args>(args)...),\
						boost::accumulators::weight = w.wi[i]*w.wi[j]);
					acc(f(-c1x*w.xi[i] + c2x, c1y*w.xi[j] + c2y, std::forward<Args>(args)...),\
						boost::accumulators::weight = w.wi[i]*w.wi[j]);
					acc(f(-c1x*w.xi[i] + c2x,-c1y*w.xi[j] + c2y, std::forward<Args>(args)...),\
						boost::accumulators::weight = w.wi[i]*w.wi[j]);
					j += 1;
				}
				i += 1;
			}
			return c1x*c1y*boost::accumulators::sum(acc);
		}

		template<unsigned order,typename RetT, typename... Args>
		RetT integrateGL(RetT (*f)(double x, Args&&... ), double from, double to, Args... args)
		{
			AccT<RetT> acc;
			GLWeights::GLWeights<order> w;
			if(! GLWeights::defined<GLWeights::GLWeights<order>>::value)
			{
				//LOG(ERROR) << "no Legendre polynomial cache for this order!";
				std::cout << "error, order not cached" << std::endl;
			}
			const double c1 = (to-from)/2.0;
			const double c2 = (to+from)/2.0;

			int i = 0;
			if(order%2)
			{
				const int lastEl =order/2;
				acc(f(c2, std::forward<Args>(args)...), boost::accumulators::weight = w.wi[0]);
				i += 1;
			}
			while(i<(order+1)/2)
			{
				acc(f( c1*w.xi[i] + c2, std::forward<Args>(args)...), boost::accumulators::weight = w.wi[i]);
				acc(f(-c1*w.xi[i] + c2, std::forward<Args>(args)...), boost::accumulators::weight = w.wi[i]);
				i += 1;
			}
			return c1*boost::accumulators::sum(acc);
		}
		template<unsigned order,typename RetT, typename... Args>
		RetT integrateRiemann(RetT (*f)(double, Args&&... ), double from_x, double to_x, Args... args)
		{
			RetT res = 0.0;
			for(double x=from_x;x<to_x;x+= (to_x-from_x)/order){
				res += f(x,std::forward<Args>(args)...);
			}
			res = res/(order);
			return res;
		}

		template<unsigned order,typename RetT, typename... Args>
		RetT integrateRiemann(RetT (*f)(double, double, Args&&... ), double from_x, double to_x, double from_y, double to_y, Args... args)
		{
			RetT res = 0.0;
			for(double y=from_y;y<to_y;y+= (to_y-from_y)/order){
			for(double x=from_x;x<to_x;x+= (to_x-from_x)/order){
				res += f(x,y,std::forward<Args>(args)...);
			}
			}
			res = res/(order*order);
			return res;
		}

};



} //end namespace utility
#endif
