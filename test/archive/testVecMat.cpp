#include <iostream>
#include <array>
#include "lettuce/core/Vec.h"
#include "lettuce/core/Mat.h"

bool testVector ();
bool testMatrix ();
void PrintResult(bool result, const std::string& label);

int main() {
    
    PrintResult(testVector(), "Vector");
    PrintResult(testMatrix(), "Matrix");

    return 0;
}


bool testVector (){

    Vec<int, 3> v1 ({1, 0, 1});
    Vec<int, 3> v2 ({1, -1, 1});
    Vec<int, 3> v3 = -v1 + 2*v2;
    v3 -= (v1.dot(v2))*v2;

    Vec<int, 3> result ({-1, 0, -1});

    return !(v3 != result || v3.abs2()!= 2);
}

bool testMatrix (){

    Vec<double, 2> vec1({1, -1});  
    Vec<double, 2> vec2({2, 1});
    Mat<double, 2, 2> mat1({vec1, vec2});

    Vec<double, 2> vec3({1, 0});  
    Vec<double, 2> vec4({0, -1});
    Mat<double, 2, 2> mat2({vec3, vec4});   

    Mat<double, 2, 2> identity_mat ({vec3, -vec4});   

    return mat1*mat1.inverse() == identity_mat;
}

void PrintResult(bool result, const std::string& label){

    if (!result){
        std::cout << "there is a problem in "<< label <<" implementation!" << std::endl;
    }else{
        std::cout << label << " test passed."<< std::endl;
    }
}
