#include <iostream>
#include <array>
#include "Vec.h"
#include "Mat.h"

void test_vector ();
void test_matrix ();

int main() {
    
    test_vector();
    test_matrix();

    return 0;
}


void test_vector (){

    Vec<int, 3> v1 ({1, 0, 1});
    Vec<int, 3> v2 ({1, -1, 1});
    Vec<int, 3> v3 = -v1 + 2*v2;
    v3 -= (v1.dot(v2))*v2;

    Vec<int, 3> result ({-1, 0, -1});

    if (v3 != result || v3.abs2()!= 2){
        std::cout << "there is a problem in Vector implementation!" << std::endl;
    }else{
        std::cout << "Vector test passed."<< std::endl;
    }
}

void test_matrix (){

    Vec<double, 2> vec1({1, -1});  
    Vec<double, 2> vec2({2, 1});
    Mat<double, 2, 2> mat1({vec1, vec2});

    Vec<double, 2> vec3({1, 0});  
    Vec<double, 2> vec4({0, -1});
    Mat<double, 2, 2> mat2({vec3, vec4});   

    Mat<double, 2, 2> identity_mat ({vec3, -vec4});   

    if (mat1*mat1.inverse() != identity_mat){
        std::cout << "there is a problem in Matrix implementation!" << std::endl;
    }else{
        std::cout << "Matrix test passed."<< std::endl;
    }
}
