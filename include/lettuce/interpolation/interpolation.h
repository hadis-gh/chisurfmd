#pragma once

#include "lettuce/utils.h"

namespace lettuce
{

template<class Interpolation, class Hold, typename... TArgs>
auto interpolate(Interpolation interpolator, const Hold& hold, TArgs... xyz)
{
	return interpolator(hold, xyz...);
}

template<class Interpolation>
auto interpolate(Interpolation interpolator, double y)
{
	return y;
}

template<class Interpolation>
auto interpolate(Interpolation interpolator, float y)
{
	return y;
}

struct InterpolateNearest
{
	template<class Hold, typename T, typename... TArgs>
	T operator()(const Hold& hold, T x, TArgs... yz)
	{
		auto idx = hold.getClosestIdx(x);
		return static_cast<T>(interpolate(*this, idx.y(), yz...));
	};
};

struct InterpolateLinear
{
	template<class Hold, typename T, typename... TArgs>
	T operator()(const Hold& hold, T x, TArgs... yz)
	{
		auto idx = hold.getClosestIdx(x);
		const T y1 = interpolate(*this, idx.y(), yz...);
		const T x1 = idx.x();
		if(lettuce::floateq(x1, x))
			return y1;

		const auto idx2 = (idx.x() < x) ? ++idx : --idx;
		if(idx == idx2)
			return y1;

		const T y2 = interpolate(*this, idx2.y(), yz...);

		const T x2 = idx2.x();

		return ( y1*(x2 - x) + y2*(x - x1) ) / (x2 - x1);
	};
};

}
