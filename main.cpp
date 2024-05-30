#include <iostream>
#include <array>
#include "Vec.h"
#include "Mat.h"

int main() {
    
    Vec<int, 3> v1 ({1, 0, 1});
    Vec<int, 3> v2 ({1, -1, 1});
    Vec<int, 3> v3 ({0, 0, 1});
    Vec<double, 3> v3d ({0, 0, 1});
    std::cout << "-------------------"<< std::endl;
    std::cout << "Vector test"<< std::endl;
    std::cout << "-------------------"<< std::endl;
    v3 += v1;
    std::cout << v3 << std::endl;
    v3 -= v1;
    std::cout << v3 << std::endl;
    v3d -= v1;
    std::cout << v3d << std::endl;
    v3 *= 5.;
    std::cout << v3 << std::endl;
    std::cout << -v1 << std::endl;
    std::cout << v1 + v2 << std::endl;
    std::cout << v1 - v2 << std::endl;
    std::cout << v1 * v2 << std::endl;
    std::cout << v1 * 2 << std::endl;
    std::cout << 2 * v1 << std::endl;
    std::cout << v1 * v2 << std::endl;
    std::cout << v1.dot(v2) << std::endl;    
    std::cout << v1.abs2() << std::endl;    
    std::cout << "-------------------"<< std::endl;
    std::cout << "Matrix test"<< std::endl;
    std::cout << "-------------------"<< std::endl;

    // Mat<int, 2, 2> mat1 = {{1, 0}, {0, 1}};




    return 0;
}