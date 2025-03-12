#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <iterator>

#include "Traits.h"

namespace lettuce
{

template<class T, class Y,
	template<class> class YContainer = std::vector, template<class> class XContainer = std::vector>
class IrregularSelector
{
	XContainer<T> m_x;
	YContainer<Y> m_y;

#if 0
	std::pair<unsigned int, unsigned int> getNeighborsIdxInner(T x) const noexcept
	{
		// linear search
		for(unsigned int a = 1; a < m_x.size(); ++a)
		{
			if(x < m_x[a])
				return {a-1, a};
		}

		// TODO binary search?
		// auto pivot = m_x.size() / 2;
	}
#endif

	public:

	using YContainerType = YContainer<Y>;
	using XContainerType = XContainer<T>;
	using XType = T;
	using YType = Y;

	const XContainerType& x() const {return m_x;}
	const YContainerType& y() const {return m_y;}
	std::size_t nPoints() const {return m_x.size();}

	class Idx
	{
		protected:
		const IrregularSelector<T,Y, YContainer, XContainer>& m_cont;
		unsigned int m_idx;

		public:

		Idx(const IrregularSelector<T, Y, YContainer, XContainer>& cont, unsigned int idx)
			: m_cont(cont), m_idx(idx)
		{}

		decltype(auto) y() const noexcept {return m_cont.m_y[m_idx];}
		T x() const noexcept {return m_cont.m_x[m_idx];}

		Idx operator++() const noexcept
		{
			const auto idx = m_idx + 1;
			return Idx(m_cont, (idx == m_cont.m_x.size()) ? m_idx : idx);
		}

		Idx operator--() const noexcept
		{
			const auto idx = (m_idx == 0) ? m_idx : m_idx - 1;
			return Idx(m_cont, idx);
		}
		
		bool atBegin() const noexcept { return m_idx == 0; }
		bool atEnd() const noexcept { return m_idx == m_cont.m_x.size()-1; }

		bool operator==(const Idx& other) const noexcept
		{
			return m_idx == other.m_idx;
		}

		bool operator!=(const Idx& other) const noexcept
		{
			return m_idx != other.m_idx;
		}
	};

	IrregularSelector(XContainer<T> x, YContainer<Y> y)
		: m_x(std::move(x)), m_y(std::move(y))
	{
		if(m_x.empty())
		{
			throw std::runtime_error("IrregularSelector cannot be constructed with empty set.");
		}
		if(m_x.size() != m_y.size())
		{
			std::ostringstream os;
			os << "IrregularSelector requires x and y arrays to be of the same size, got "
				"x: " << m_x.size() << " y: " << m_y.size();
			throw std::runtime_error(std::move(os).str());
		}

		auto i = m_x.begin();
		auto prev = *i;
		for(++i; i != m_x.end(); ++i)
		{
			if(prev >= *i)
			{
				std::ostringstream os;
				os << "IrregularSelector selector requires strictly ascending ordinates. Violated at "
					<< std::distance(m_x.begin(), i) << " by !( " << prev << " < " << *i << ").";
				throw std::runtime_error(std::move(os).str());
			}
			prev = *i;
		}
	}

	Idx getClosestIdx(T x) const noexcept
	{
#if 0
		const auto it = std::find_if(m_x.begin(), m_x.end(), [x](mx){return mx > x};);
#else
		const auto it = std::upper_bound(m_x.begin(), m_x.end(), x);
#endif
		unsigned int idx = m_x.size()-1;
		if (it != m_x.end())
		{
			idx = std::distance(m_x.begin(), it);
			if(idx > 0)
				if(m_x[idx] - x > x - m_x[idx-1])
					idx -= 1;
		}
		return Idx(*this, idx);
	}

#if 0
	std::pair<unsigned int, unsigned int> getNeighborsIdx(T x) const noexcept
	{
		if (x <= m_x.front())
			return {0, 1};
		if (x >= m_x.back())
			return {m_x.size() - 2, m_x.size() - 1};

		return getNeighborsIdxInner(x);
	}
#endif
};

namespace traits
{
	template<class T, class Y, template<class> class YContainer, template<class> class XContainer>
	struct YType<IrregularSelector<T, Y, XContainer, YContainer>>
	{
		using type = Y;
	};
}

template<class T, class Y, template<class> class YContainer = std::vector, template<class> class XContainer = std::vector>
auto makeIrregularSelector(XContainer<T> x, YContainer<Y> y)
{
	return IrregularSelector<T, Y>(x, y);
}

}
