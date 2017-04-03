#ifndef MEMOIZATION_H_
#define MEMOIZATION_H_

#include <boost/type_traits.hpp>	// this seems a bit better, than std::type_traits
#include <boost/unordered_map.hpp>	// std::unordered_map does not support incomplete types
// tests suggest unordered_map for static, map for dynamic
#include <cstddef>			// nullptr_t
#include <functional>	    // less,hash
#include <utility>			// std::move, std::forward
#include <algorithm>		// std::find
#include <iterator>			// std::cbegin, cend

/*!
 *	Goals and features for this memoization library:
 *
 *	= template parameters: function pointer, return type, table size
 *
 *	= dynamic
 *		- table size defines max size. either override on full or call function 
 *
 *	= static (not implemented yet)
 *		- build a static lookup table (hash)
 *		- directly hash arguments to next table entry (interpolation)
 *		- complete vs incomplete:
 *			- for incomplete table include fallback to direct function call
 *
 *	= specialized (not implemented yet)
 *		- Matsubara Freq
 *		- exp(ix)	
 *
 *	= possible issues:
 *		- find a way to pass required alignments (specify or track access pattern)
 *		- this accepts almost everything (function, functor, lambda), there really should be some checks (see ref -1- for more)
 *		- there seems to be no obvious way in C++14 to check whether a given function is pure
 *		- if std::function or boost:function becomes more usable, much of this would simplify a lot (FPointType vanishes)
 * 
 * 
 *	this takes inspiration from (especially -1- and -5-):
 *		-1- Jim Porter, memo: https://github.com/jimporter/memo
 *		-2- http://www.bigoh.co.uk/c++/2016/02/06/memoization.html
 *		-3- http://cpptruths.blogspot.de/2012/01/general-purpose-automatic-memoization.html
 *		-4- https://projects.giacomodrago.com/c++memo
 *		-5- Scott Meyers: Effective Modern C++
 */

namespace Memoize
{
using std::cbegin;
using std::cend;
using std::find;

// taken from source -2-
// detail.function_signature is set by template matching
namespace detail {
	// This gets the base function type from pretty much anything it can:
	// C function types, member function types, monomorphic function objects.
	template<typename T, typename Enable = void>
	struct function_signature;

	// match functor (overloaded operator() for class).  
	template<typename T>
	struct function_signature<T, typename std::enable_if<std::is_class<T>::value>::type> : 
		function_signature<decltype(&T::operator())> {};

	// member function of T
	template<typename T, typename Ret, typename... Args>
	struct function_signature<Ret (T::*)(Args...)> {
		using type = Ret(Args...);
	};

	template<typename T, typename Ret, typename... Args>
	struct function_signature<Ret (T::*)(Args...) const> {
		using type = Ret(Args...);
	};

	// standard functions (f pointers and refs)
	template<typename Ret, typename... Args>
	struct function_signature<Ret(Args...)> {
		using type = Ret(Args...);
	};

	template<typename Ret, typename... Args>
	struct function_signature<Ret (&)(Args...)> {
		using type = Ret(Args...);
	};

	template<typename Ret, typename... Args>
	struct function_signature<Ret (*)(Args...)> {
		using type = Ret(Args...);
	};
}


template<typename FPointType, typename Function>
class Memoizer;
/*!	 @brief	FPointType matches the function type, templates for almost all possible types are taken from source -2- 
 */
template <typename FPointType, typename ReturnType, typename... ArgTypes>
class Memoizer<FPointType, ReturnType(ArgTypes...)>
{

    private:
		// only get bare types, see -5- item 9
		using ArgTupleType  = std::tuple<typename boost::remove_reference<ArgTypes>::type...>;
		//using MapType	    = std::unordered_map<FPointType>;
		using MapType	    = std::map<ArgTupleType, ReturnType, std::less<void>>;
		using HMapType		= boost::unordered::unordered_map<ArgTupleType, ReturnType>;
		using ReturnRef	    = const ReturnType&;


		MapType	    _mem;   // Storage for results
		FPointType  _f;	    // Storage for function pointer
		//TODO: placeholder 
		std::nullptr_t _hf;	// Storage for hash function pointer (only used for hash map storage)
		// -1- does a check here in a private sub class whether the given type can actually be memoized
		// this is omitted for now

    public:
		Memoizer(const FPointType &f, HMapType hmap, const size_t N = 0): _f(f), _hf(hf) {}	
		Memoizer(const FPointType &f): Memoizer(f,0) {}	

		/*!	@brief	encapsulates memoization. For versions which can't be found in 
		 *			storage (overload for euqal possible), generate new entry (if size < N or N == 0)
		 */
		auto operator()(ArgTypes... args) -> ReturnRef
		{
			// faster version of ... = make_tuple(args...); see -5- items 23, 28
			// prefere const iterators over normal ones this only works for C++14, see -5-, item 42
			// BUG: intel compiler has not fully implemented C++14 (std::cbegin cen still missing from algorithm)
			const auto argsTuple = std::forward_as_tuple(std::forward<ArgTypes>(args)...);
			auto memPtr = _mem.find(argsTuple);
			//std::find(_mem.begin(),_mem.end(),argsTuple);
			if(memPtr != std::cend(_mem)){						// found args tuple in _mem
				return memPtr->second;						// return stored val
			}else{
				auto const result = _f(args...);
				// return reference to inserted element, use emplace to avoid unecessary constructor calls
				// during insertion with moved (->RValue) arguments. see -5- items 23 and 42
				return _mem.emplace(std::move(argsTuple),std::move(result)).first->second;
			}
		}
};

/*!
 *	basis teplate for memeoize(). Used after matching Memoize::detail
 *	this uses universal references (see -5- item 24)
 */
template<typename FPointType, typename Function >
inline auto memoize(FPointType&& fp) {
	return Memoizer<FPointType, Function>(std::forward<FPointType>(fp));
}

/*!	@brief iconstructs and returns function handle with memoizing backend
 *
 *	@param fp RValue reference to function 
 *
 *	@return	function pointer with the same outputs as the input function, but internal memoizing
 */
template<typename FPointType>
inline auto memoize(FPointType&& fp) {
	return Memoizer<FPointType, typename detail::function_signature<FPointType>::type>(	
			std::forward<FPointType>(fp) 
		);
}

/*!	@brief	
 *
 *
 */
template<typename FPointType, typename HashFPointType>
inline auto memoize(FPointType&& fp, HashFPointType&& hfp) {
	return Memoizer<FPointType, typename detail::function_signature<FPointType>::type>(	
			std::forward<FPointType>(fp) 
		);
}

}
#endif
