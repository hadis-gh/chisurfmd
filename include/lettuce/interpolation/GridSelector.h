#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

#include "lettuce/utils.h"

#include "Traits.h"

namespace lettuce
{

template<typename T>
auto getCheckGrid(const std::vector<T>& x)
{
	if(x.size() < 2)
		throw std::runtime_error("Cannot deduce grid from <2 points.");

	const auto delta = x[1] - x[0];
	for(std::size_t a = 2; a < x.size(); ++a)
	{
		const auto tdelta = x[a] - x[a-1];
		if(!floateq(delta, tdelta))
		{
			std::ostringstream os;
			os << "Data points not regular. Expected " << delta
				<< ", but found " << tdelta << " at " << a << " - " << a-1 << ".";
			throw std::runtime_error(std::move(os).str());
		}
	}

	return delta;
}

template<class T, class Y, template<class> class YContainer = std::vector, bool Periodic = false>
class GridSelector
{
	T m_min, m_step;
	unsigned int m_last;
	YContainer<Y> m_y;

	public:

	using YContainerType = YContainer<Y>;
	using XType = T;
	using YType = Y;
	static constexpr bool isPeriodic() {return Periodic;}

	const YContainerType& y() const {return m_y;}
	T xmin() const {return m_min;}
	T xstep() const {return m_step;}
	T xlast() const {return m_last;}
	std::size_t nPoints() const {return m_last+1;}

	class Idx
	{
		protected:
		const GridSelector<T,Y, YContainer, Periodic>& m_cont;
		unsigned int m_idx;

		public:

		Idx(const GridSelector<T, Y, YContainer, Periodic>& cont, unsigned int idx) noexcept
			: m_cont(cont), m_idx(idx)
		{}

		decltype(auto) y() const noexcept {return m_cont.m_y[m_idx];}
		T x() const noexcept {return m_cont.m_min + m_cont.m_step*m_idx;}

		Idx operator++() const noexcept
		{
			unsigned int idx;
			if constexpr (Periodic)
				idx = (m_idx == m_cont.m_last) ? 0 : m_idx + 1;
			else
				idx = (m_idx == m_cont.m_last) ? m_idx : m_idx + 1;
			return Idx(m_cont, idx);
		}

		Idx operator--() const noexcept
		{
			unsigned int idx;
			if constexpr (Periodic)
				idx = (m_idx == 0) ? m_cont.m_last : m_idx - 1;
			else
				idx = (m_idx == 0) ? m_idx : m_idx - 1;
			return Idx(m_cont, idx);
		}
		
		bool atBegin() const noexcept { return m_idx == 0; }
		bool atEnd() const noexcept { return m_idx == m_cont.m_last; }

		bool operator==(const Idx& other) const noexcept
		{
			return m_idx == other.m_idx;
		}

		bool operator!=(const Idx& other) const noexcept
		{
			return m_idx != other.m_idx;
		}
	};

	GridSelector(T min, T step, YContainer<Y> y)
		: m_min(min), m_step(step), m_y(std::move(y))
	{
		m_last = m_y.size()-1;
	}

	GridSelector(T min, T step, unsigned int last, YContainer<Y> y)
		: m_min(min), m_step(step), m_last(last), m_y(std::move(y))
	{
	}

	template<class XContainer>
	GridSelector(const XContainer& x, YContainer<Y> y)
		: m_y(std::move(y))
	{
		m_step = getCheckGrid(x);
		m_min = x[0];
		m_last = x.size()-1;
	}

	Idx getClosestIdx(T x) const noexcept
	{
		if constexpr (Periodic)
		{
			const auto xsup = m_step*(m_last + 1);
			x = std::fmod(x - m_min, xsup);
			if(x < 0)
				x += xsup;
			x += m_min;
			unsigned int idx = std::round((x - m_min)/m_step);
			if (idx > m_last)
				idx = 0;
			return Idx(*this, idx);
		}
		else
		{
			unsigned int idx = 0u;
			if (x > m_min)
				idx = std::round((x - m_min)/m_step);
			return Idx(*this, std::min(idx, m_last));
		}
	}

#if 0
	std::pair<unsigned int, unsigned int> getNeighborsIdx(T x) const noexcept
	{
		if (x < m_min)
			return {0, 1};
		const T idx = (x - m_min)/m_step;
		if(idx >= m_last)
			return {m_last-1, m_last};
		return {static_cast<unsigned int>(std::floor(idx)), static_cast<unsigned int>(std::ceil(idx))};
	}
#endif
};

namespace traits
{
	template<class T, class Y, template<class> class YContainer>
	struct YType<GridSelector<T, Y, YContainer>>
	{
		using type = Y;
	};
}


template<class T, class Y, template<class> class YContainer = std::vector>
auto makeGridSelector(T min, T step, YContainer<Y> y)
{
	return GridSelector<T, Y>(min, step, y);
}

template<class T, class Y, template<class> class YContainer = std::vector>
auto makeGridSelector(T min, T step, unsigned int last, YContainer<Y> y)
{
	return GridSelector<T, Y>(min, step, last, y);
}

}
