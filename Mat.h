#pragma once

#include <array>
#include <iostream>
#include "Vec.h"

template<typename T, int M = 2, int N = 2>
class Mat {
    public:
    std::array <Vec<T, N>, M> elements;

    Mat() = default;
    Mat(std::array<Vec<T, N>, M> init) : elements(init) {}

    T& operator()(int i, int j);
    const T& operator()(int i, int j) const;

    Vec<T, N>& operator()(int i);
    const Vec<T, N>& operator()(int i) const;

    Mat<T, M, N> transpose() const;
    
    template <int MM, int NN>
    Mat<T, M, N> operator* (const Mat<T, MM, NN> &other) const;
    Mat<T, M, N> &operator* (const T &scaler);

    template <int MM, int NN>
    Mat<T, M, N> inner_product (const Mat<T, MM, NN> &other);

    template <int MM, int NN>
    Mat<T, M, N> outer_product (Mat<T, MM, NN> &other) const;

    T get_determinant() const;

    Mat<T, M, N> inverse() const;

    template <int MM, int NN>
    bool operator==(const Mat<T, MM, NN> &other) const;

    template <int MM, int NN>
    bool operator!=(const Mat<T, MM, NN> &other) const;   
};

template<typename T, int M = 2, int N = 2>
std::ostream& operator<<(std::ostream& COUT, const Mat<T, M, N>& ref);

#include "Mat.inl"