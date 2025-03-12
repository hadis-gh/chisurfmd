#pragma once

#include <type_traits>

namespace lettuce
{

namespace traits
{
	template<class SelectorType, typename SFINAE = void>
	struct YType;

	template<class T>
	struct YType<T>
	{
		using type = T;
	};

	template<class SelectorType, typename SFINAE = void>
	struct SelectorStride
	{
		static constexpr std::size_t stride(const SelectorType&)
		{
			return 1u;
		}
	};
}

template<class T>
using YType = traits::YType<T>::type;

template<class SelectorType>
std::size_t selectorStride(const SelectorType& s)
{
	return traits::SelectorStride<SelectorType>::stride(s);
}

namespace detail
{
	template<class Hold, std::size_t C = 1>
	struct HoldDims
	{
		static constexpr auto value = HoldDims<typename Hold::YType, C+1>::value;
	};

	template<class Hold, std::size_t C>
		requires std::is_floating_point_v<typename Hold::YType>
	struct HoldDims<Hold, C>
	{
		static constexpr std::size_t value = C;
	};
}

template<class Hold>
constexpr std::size_t HoldDims_v = detail::HoldDims<Hold>::value;

}
