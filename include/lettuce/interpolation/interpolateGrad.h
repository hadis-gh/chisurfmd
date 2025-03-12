#pragma once

#include <array>
#include <tuple>
#include <type_traits>
#include <iostream>

#include "lettuce/utils.h"

#include "interpolation.h"

namespace lettuce
{

namespace detail
{

template<class Hold, int C, typename... XYZ>
struct InterpolateGradLinear;

template<class Hold, int C, typename X, typename... YZ>
struct InterpolateGradLinear<Hold, C, X, YZ...>
{
	static auto get(const Hold& hold, X x, YZ... yz)
	{
		auto idx = hold.getClosestIdx(x);
		auto grad = InterpolateGradLinear<std::decay_t<decltype(idx.y())>, C+1, YZ...>::
			get(idx.y(), std::forward<YZ>(yz)...);
		
		const auto y = grad[C];
		using T = decltype(y);
		auto ng = [&](auto nidx){
			const auto ny = interpolate(InterpolateNearest(), nidx.y(), std::forward<YZ>(yz)...);
			// std::cout << ny << '\t';
			if(floateq(nidx.x(), idx.x()))
				return std::make_tuple(0., double(idx.x()));
			else
				return std::make_tuple(double((ny - y) / (nidx.x() - idx.x())), double( nidx.x() + idx.x() ) / 2.);
		};

		const auto [y1, x1] = ng(--idx);
		const auto [y2, x2] = ng(++idx);
		// std::cout << "grad[" << C << "] " << y1 << '\t' << y << '\t' << y2
		// 	<< '\t' << x1 << '\t' << x << '\t' << x2 << '\n';

		if(floateq(x1, x2) || (x < x1) || (x > x2))
			grad[C] = 0.;
		else
			grad[C] = ( y1*(x2 - x) + y2*(x - x1) ) / (x2 - x1);

		return grad;
	};
};

template<class T, int C>
	requires std::is_floating_point_v<T>
struct InterpolateGradLinear<T, C>
{
	static std::array<T, C> get(const T& v)
	{
		std::array<T, C> g;
		g.fill(v);
		return g;
	};
};

}

template<class Hold, typename... XYZ>
auto interpolateGradLinear(const Hold& hold, XYZ... xyz)
{
	return detail::InterpolateGradLinear<Hold, 0, XYZ...>::get(hold, std::forward<XYZ>(xyz)...);
}

}
