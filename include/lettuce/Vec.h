#pragma once

#include <iostream>
#include <array>

template<typename T , int D = 2> 
class Vec {
    std::array<T, D> elements;

public:
    T x = elements[0];
    T y = elements[1];
    T z = elements[2];
    
    // Vec() = default;
    Vec(T v=0) {elements.fill(v);}
    Vec(std::array<T, D> init) : elements(init) {}
               
    T& operator[](int index);
	const T& operator[] (int index) const;

    Vec<T, D> operator+(const Vec<T, D>& other) const;
    Vec<T, D> operator-(const Vec<T, D>& other) const;
    Vec<T, D> operator-() const;
    
    T operator* (const Vec<T, D>& other) const;
    Vec<T, D> operator*(const T& number) const;

    template<typename TT>
    Vec<T, D>& operator+=(const Vec<TT, D>& other); 
    template<typename TT>    
    Vec<T, D>& operator-=(const Vec<TT, D>& other);
    Vec<T, D>& operator *= (const T& number);

    template<typename TT>    
    bool operator==(const Vec<TT, D>& other)const;

    template<typename TT>    
    bool operator!=(const Vec<TT, D>& other){ return !(*this == other);}

    bool operator>=(const T& number) const;
    bool operator<=(const T& number)const;

    bool operator>(const T& number)const;
    bool operator<(const T& number)const;

    Vec<T, D> operator- (const T& number)const;
    Vec<T, D> operator+ (const T& number)const;

    T sum()const;
    T dot(const Vec<T, D>& other)const;
    T abs2()const;
};

template<typename T, int D>
Vec<T, D> operator*(const T& number, const Vec<T, D>& vec);

template<typename T, int D>
std::ostream& operator<<(std::ostream& COUT, const Vec<T, D>& ref);

#include "Vec.inl"