#pragma once

#include <sstream>
#include <vector>
#include <stdexcept>
#include <utility>
#include <cmath>
#include <algorithm>
#include <deque>
#include <span>
#include <list>

#include "lettuce/utils.h"

#include "GridSelector.h"
#include "IrregularSelector.h"
#include "UnifyingYContainer.h"

namespace lettuce
{

template<typename T>
auto getXRanges(std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
{
	std::span<T> rspan(table[c].begin() + beginr, table[c].begin() + endr);

	std::vector<T> x;
	std::vector<std::size_t> ranges;
	{
		x.reserve(rspan.size());
		ranges.reserve(rspan.size());
		ranges.push_back(0);
		auto prev = rspan.front();
		x.push_back(prev);
		for(std::size_t r = 1; r < rspan.size(); ++r)
		{
			if(prev != rspan[r])
			{
				ranges.push_back(r);
				prev = rspan[r];
				x.push_back(prev);
			}
		}
		x.shrink_to_fit();
		ranges.push_back(rspan.size());
		ranges.shrink_to_fit();
	}

	// std::cout << std::span(x) << x.size() << '\n';
	// std::cout << std::span(ranges) << ranges.size() << '\n';

	return std::make_tuple(std::move(rspan), std::move(x), std::move(ranges));
}

namespace fillHoldDetails
{
	template<class Hold, typename SFINAE = void>
	struct FillHoldInner;

template<class T>
void deepSortRange(std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
{
	std::vector<std::size_t> idxlist(endr-beginr);
	for(std::size_t i = 0; i < idxlist.size(); ++i)
		idxlist[i] = i;
	std::span<T> rspan(table[c].begin() + beginr, table[c].begin() + endr);
	std::stable_sort(idxlist.begin(), idxlist.end(), [rspan](auto a, auto b){return rspan[a] < rspan[b];});

	{
		std::vector<T> x(idxlist.size());
		for(int cc = c; cc < table.size(); ++cc)
		{
			for(int a = 0; a < idxlist.size(); ++a)
				x[a] = table[cc][beginr + idxlist[a]];
			std::copy(x.cbegin(), x.cend(), table[cc].begin()+beginr);
		}
	}
}

template<class Hold>
struct FillHold
{
	template<class T = float>
	static Hold make(std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		deepSortRange(table, c, beginr, endr);

		return FillHoldInner<Hold>::make(table, c, beginr, endr);
	}
};

template<class Hold, typename SFINAE = void>
struct FillHoldUnified
{
	template<class T = float>
	static Hold unify(std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr
			, const std::vector<std::size_t>& ranges)
	{
		using InnerYType = typename Hold::YType;
		std::vector<InnerYType> innerY;
		std::vector<T> unifiedX;
		for(int r = 1; r < ranges.size(); ++r)
		{
			// std::cout << r << "\tr/" << ranges.size() << '\t' <<ranges[r-1] << '\t' << ranges[r] << std::endl;
			deepSortRange(table, c, ranges[r-1], ranges[r]);
			const auto beginrr = beginr+ranges[r-1];
			auto [rspan, x, innerRanges] = getXRanges(table, c, beginrr, beginr+ranges[r]);
			if(r == 1)
			{
				unifiedX = std::move(x);
				innerY.reserve((ranges.size()-1) * unifiedX.size());
			}
			else
			{
				for(std::size_t a = 0; a < x.size(); ++a)
				{
					if(x[a] != unifiedX[a])
					{
						std::ostringstream os;
						os << "Cannot unify dimension " << c << " found x-mismatch at index "
							<< a << " : " << unifiedX[a] << " != "  << x[a] << " .";
						throw std::runtime_error(std::move(os).str());
					}
				}
			}

			for(int ir = 1; ir < innerRanges.size(); ++ir)
			{
				// std::cout << ir << "\t/" << innerRanges.size() << '\t' << innerRanges[ir-1] << '\t' << innerRanges[ir] << std::endl;
				innerY.emplace_back(FillHold<InnerYType>::make(table, c+1
							, beginrr + innerRanges[ir-1], beginrr + innerRanges[ir]));
			}
		}

		return {std::move(unifiedX), std::move(innerY)};
	}
};

template<class Hold>
	requires (std::is_same_v<
			typename Hold::YContainerType, UnifyingYContainer<typename Hold::YType>>)
struct FillHoldUnified<Hold>
{
	template<class T = float>
	static Hold make(std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		deepSortRange(table, c, beginr, endr);

		return FillHoldInner<Hold>::make(table, c, beginr, endr);
	}

	template<class T = float>
	static Hold unify(std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr
			, const std::vector<std::size_t>& ranges)
	{
		using InnerYType = typename Hold::YType;
		std::vector<T> unifiedX;
		std::list<std::vector<std::size_t>> allInnerRangesList;
		for(int r = 1; r < ranges.size(); ++r)
		{
			deepSortRange(table, c, ranges[r-1], ranges[r]);
			const auto beginrr = beginr+ranges[r-1];
			auto [rspan, x, innerRanges] = getXRanges(table, c, beginrr, beginr+ranges[r]);
			if(r == 1)
			{
				unifiedX = std::move(x);
			}
			else
			{
				for(std::size_t a = 0; a < x.size(); ++a)
				{
					if(x[a] != unifiedX[a])
					{
						std::ostringstream os;
						os << "Cannot unify dimension " << c << " found x-mismatch at index "
							<< a << " : " << unifiedX[a] << " != "  << x[a] << " .";
						throw std::runtime_error(std::move(os).str());
					}
				}
			}

			for(auto& ir : innerRanges)
				ir += beginrr;

			// std::cout << "push back range ";
			// if(!allInnerRangesList.empty())
			// 	std::cout << allInnerRangesList.back().back() << '\t';
			// std::cout << innerRanges.front() << std::endl;
			allInnerRangesList.push_back(innerRanges);
		}

		const std::size_t nAllInnerRanges = [&allInnerRangesList](){
			std::size_t nAllInnerRanges = 0;
			for(const auto& v: allInnerRangesList)
				nAllInnerRanges += v.size() - 1;
			return nAllInnerRanges + 1;
			}();

		std::vector<std::size_t> allInnerRanges(nAllInnerRanges);
		auto beginIr = allInnerRanges.begin();
		for(const auto& v: allInnerRangesList)
		{
			std::copy(v.begin(), v.end(), beginIr);
			beginIr += v.size() - 1;
		}
		// std::cout << "allInnerRanges.size()= " << allInnerRanges.size() << "\t" << allInnerRanges[0] << "\t" << allInnerRanges[1] << std::endl;

		return {std::move(unifiedX),
				FillHoldUnified<InnerYType>::unify(table, c+1
					, beginr, endr, allInnerRanges)
			};
	}
};

// template<typename T, class Y,
// 	template<class> class YContainer, template<class> class XContainer>
template<class Hold>
struct FillHoldInner<Hold>
{
	template<typename T>
	static Hold make(
			std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		auto [rspan, x, ranges] = getXRanges(table, c, beginr, endr);

		typename Hold::YContainerType y;
		y.reserve(x.size());

		for(int r = 1; r < ranges.size(); ++r)
		{
			y.emplace_back(FillHold<typename Hold::YType>::make(table
						, c+1, beginr+ranges[r-1], beginr+ranges[r]));
		}

		return {std::move(x), std::move(y)};
	}
};

// template<template<class...> Selector, typename T, class Y>
template<class Hold>
	requires std::is_same<typename Hold::YContainerType, UnifyingYContainer<typename Hold::YType>>::value
struct FillHoldInner<Hold>
{
	// using T = typename Hold::XType;
	using Y = typename Hold::YType;

	template<typename T>
	static Hold make(
			std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		auto [rspan, x, ranges] = getXRanges(table, c, beginr, endr);

		UnifyingYContainer<Y> y(FillHoldUnified<Y>::unify(table, c+1, beginr, endr, ranges));

		return {std::move(x), std::move(y)};
	}
};

template<typename T>
struct FillHoldInner<GridSelector<T, T, std::span>>
{
	static GridSelector<T, T, std::span> make(
			std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		// TODO: not all of this is needed for this case, x = rspan should do
		auto [rspan, x, ranges] = getXRanges(table, c, beginr, endr);

		const auto delta = getCheckGrid(x);

		return {0, delta,
			std::span<T>(table[c+1].begin() + beginr, table[c+1].begin() + endr)};
	}
};

/*! Reduce input points to main grid, defined by first pair. */
template<typename T>
struct FillHoldInner<GridSelector<T, T, std::vector>>
{
	static GridSelector<T, T, std::vector> make(
			std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		std::span<T> rspan(table[c].begin() + beginr, table[c].begin() + endr);

		if(rspan.size() < 2)
			throw std::runtime_error("Cannot deduce grid from <2 points.");

		const auto delta = [&](){
				const double d = rspan[1] - (double)rspan[0];
				const double range = rspan.back() - (double)rspan.front();
				return range / std::round(range/d);
			}();
		const std::size_t steps = std::round((rspan.back() - (double)rspan.front())/delta) + 1;
		std::span<T> yspan(table[c+1].begin() + beginr, table[c+1].begin() + endr);
		std::vector<T> y(steps);
		y[0] = yspan[0];
		y[1] = yspan[1];

		for(std::size_t a = 2, ay = 2; a < rspan.size() && ay < y.size(); ++a)
		{
			if(floateq((double)rspan[a], a*delta + rspan[0]))
			{
				y[ay] = yspan[a];
				++ay;
			}
		}

		return {0, (float)delta, std::move(y)};
	}
};

template<typename T>
struct FillHoldInner<IrregularSelector<T, T, std::span, std::span>>
{
	static IrregularSelector<T, T, std::span, std::span> make(
			std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		return {std::span<T>(table[c].begin() + beginr, table[c].begin() + endr),
			std::span<T>(table[c+1].begin() + beginr, table[c+1].begin() + endr)};
	}
};

template<typename T>
struct FillHoldInner<IrregularSelector<T, T, std::span, std::vector>>
{
	static IrregularSelector<T, T, std::span, std::vector> make(
			std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
	{
		std::vector<T> x(endr - beginr);
		std::copy(table[c].begin() + beginr, table[c].begin() + endr, x.begin());
		return {std::move(x),
			std::span<T>(table[c+1].begin() + beginr, table[c+1].begin() + endr)};
	}
};

} // end fillHoldDetails

template<class Hold, class T = float>
Hold fillHold(std::vector<std::vector<T>>& table, int c, std::size_t beginr, std::size_t endr)
{
	return fillHoldDetails::FillHold<Hold>::make(table, 0, 0, table.front().size());
}

} // end lettuce
