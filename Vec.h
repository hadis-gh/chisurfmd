#pragma once

#include <iostream>
#include <array>

template<typename T , int D = 2> 
class Vec {
    std::array<T, D> elements;

public:
    //constructors
    Vec() = default;
    Vec(std::array<T, D> init) : elements(init) {}
               
    T& operator[](int index);
	const T& operator[] (int index) const;

    Vec<T, D> operator+(const Vec<T, D>& other) const;
    Vec<T, D> operator-(const Vec<T, D>& other) const;
    Vec<T, D> operator-() const;
    T operator* (const Vec<T, D>& other) const;
    Vec<T, D> operator*(const T& number) const;

    Vec<T, D>& operator+=(const Vec<T, D>& other); 
    Vec<T, D>& operator-=(const Vec<T, D>& other);
    Vec<T, D>& operator *= (const T& number);
};

template<typename T, int D>
Vec<T, D> operator*(const T& number, const Vec<T, D>& vec);

template<typename T, int D>
std::ostream& operator<<(std::ostream& COUT, const Vec<T, D>& ref);

#include "Vec.inl"

using VecF2 = Vec<float, 2>;
using VecD2 = Vec<double, 2>;