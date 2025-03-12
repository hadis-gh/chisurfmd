#pragma once

namespace lettuce
{

template<typename T>
concept Selector = requires (T t, T::Idx idx) {
	typename T::Idx;
	typename T::XType;
};

template<typename T>
concept CUnifyingYContainer = requires (T t, typename T::Idx idx) {
	typename T::Idx;
	typename T::SelectorWrapper;
};

}
