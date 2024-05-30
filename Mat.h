#pragma once

#include <array>
#include <iostream>
#include "Vec.h"

// using Vec<T, D> = vec;
template<typename T, int N = 2, int M = 2>
class Mat {
    std::array <Vec<T, N>, M> elements;

public:
    //constructors
    Mat() = default;
    Mat(std::array<Vec<T, N>, M> init) : elements(init) {}
    //Accesor
    T& operator()(int i, int j) {
        return elements[i][j];
    }

    const T& operator()(int i, int j) const {
        return elements[i][j];
    }

    Mat<Vec<T, N>, M> &transpose (){
        Mat<Vec<T, N>, M> T_mat;        //should I use for loop for initialization?
        for (int i=0; i<M; i++){
            for (int j=0; j<N; j++){
                T_mat.elements[j][i] = elements[i][j]; 
            }
        }
        return T_mat;
    }
    
    template <int MM, int NN>
    Mat<Vec<T, NN>, M> operator* (const Mat<Vec<T, NN>, MM> &other) const {
        Mat<Vec<T, NN>, M> result;
        if (N != MM){
            std::cout << "Not allowded multiplication! " << N << " does not match " << MM << "!" << std::endl;
            //maybe also other action?!
        } else {
            // Initialization
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < NN; j++) {
                    result[i][j] = 0;
                }
            }

            for (int i=0; i<M; i++){
                for (int k=0; k<NN; k++){
                    for (int j=0; j<N; j++){
                        result[i][j] += elements[i][k] * other[k][j];
                        //result[i,j] += elements[i,k] * elements[k,j];   //[i,j] should work but I put [i][j] just in case
                    }
                }
            }
        }
        return result;
    }

    template <int MM, int NN>
    Mat<Vec<T, NN>, N> inner_product (const Mat<Vec<T, NN>, MM> &other){
        return ((*this).transpose()) * other;
    }

    template <int MM, int NN>
    Mat<Vec<T, MM>, M> outer_product (const Mat<Vec<T, NN>, MM> &other){
        return (*this)*(other.transpose());
    }
    // for 2D Matrix
    T get_determinant()const{
        T determinant = elements[0][0] *elements[1][1] - elements[0][1] *elements[1][0];
        return determinant;
    }

    Mat<Vec<T, N>, M> &inverse (){

        Mat<Vec<T, N>, M> inverse_mat;
        // manually calculated because we are suppose to work with 2D...
        inverse_mat [0][0] = elements[1][1];
        inverse_mat [0][1] = - elements[0][1];
        inverse_mat [1][0] = - elements[1][0];
        inverse_mat [1][1] = elements[0][0];

        for (int i=0; i<M; i++){
            for (int j=0; j<N; j++){
                (inverse_mat.get_determinant)*elements[i][j];
            }
        }
        return inverse_mat;        
    }
};