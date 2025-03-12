#pragma once

#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "Selector.h"
#include "Traits.h"

namespace lettuce
{

template<Selector Y>
class UnifyingYContainer;

namespace unifyingYContainerDetails
{

template<class Container, typename SFINAE = void>
class IndexYContainer;

template<template<class>class Container, typename Y>
	requires std::is_same<Container<Y>, std::vector<Y>>::value || std::is_same<Container<Y>, std::span<Y>>::value
struct IndexYContainer<Container<Y>>
{
	static decltype(auto) get(const Container<Y>& c, std::size_t idx, std::size_t offset)
	{
		return c[idx + offset];
	}
};

template<typename Y>
struct IndexYContainer<UnifyingYContainer<Y>>
{
	static decltype(auto) get(const UnifyingYContainer<Y>& c, std::size_t idx, std::size_t offset)
	{
		return c[idx, offset];
	}
};

template<class Container>
decltype(auto) indexYContainer(const Container& c, std::size_t idx, std::size_t offset)
{
	return IndexYContainer<Container>::get(c, idx, offset);
}

}

template<Selector Y>
class UnifyingYContainer
{
	Y m_y;
	std::size_t m_stride;

	public:

	UnifyingYContainer(Y y)
		: m_y(std::move(y))
	{
		m_stride = selectorStride(m_y) * m_y.nPoints();
	}

	const Y& y() const {return m_y;}
	std::size_t stride() const {return m_stride;}

	class Idx : public Y::Idx
	{
		using Y::Idx::m_idx;
		using Y::Idx::m_cont;
		std::size_t m_offset;

		public:
		decltype(auto) y() const noexcept
		{
			return unifyingYContainerDetails::indexYContainer(m_cont.y(), m_idx, m_offset);
		}

		Idx(const Y& y, unsigned int idx, std::size_t offset)
			: Y::Idx(y, idx), m_offset(offset)
		{}

		Idx(const Y::Idx& idx, std::size_t offset)
			: Y::Idx(idx), m_offset(offset)
		{}

		Idx operator++() const noexcept
		{
			return {++(static_cast<const Y::Idx&>(*this)), m_offset};
		}

		Idx operator--() const noexcept
		{
			return {--(static_cast<const Y::Idx&>(*this)), m_offset};
		}
	};

	class SelectorWrapper
	{
		const Y& m_y;
		std::size_t m_offset;

		public:

		SelectorWrapper(const Y& y, std::size_t idx, std::size_t stride, std::size_t offset = 0u)
			: m_y(y), m_offset(idx*stride + offset)
		{}

		const Y& operator*() const {return m_y;}

		// forward
		auto& y() const {return m_y.y();}

		// wrap
		Idx getClosestIdx(typename Y::XType x) const noexcept
		{
			return {m_y.getClosestIdx(x), m_offset};
		}
	};

	SelectorWrapper operator[](std::size_t idx) const
	{
		return SelectorWrapper(m_y, idx, m_stride);
	}

	SelectorWrapper operator[](std::size_t idx, std::size_t offset) const
	{
		return SelectorWrapper(m_y, idx, m_stride, offset);
	}
};

namespace traits
{

	template<class SelectorType>
		requires std::is_same<UnifyingYContainer<typename SelectorType::YType>, typename SelectorType::YContainerType>::value
	struct SelectorStride<SelectorType>
	{
		static std::size_t stride(const SelectorType& s)
		{
			return s.y().stride();
		}
	};

}

}
