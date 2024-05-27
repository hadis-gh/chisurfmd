#include <iostream>
#include <array>
#include "Vec.h"



int main() {
    
    Vec<int, 3> v1 ({1, 0, 1});
    Vec<int, 3> v2 ({1, -1, 1});
    Vec<int, 3> v3 ({0, 0, 1});
    v3 += v1;
    std::cout << v3 << std::endl;
    v3 -= v1;
    std::cout << v3 << std::endl;
    v3 *= 5;
    std::cout << v3 << std::endl;
    std::cout << "-------------------"<< std::endl;
    std::cout << -v1 << std::endl;
    std::cout << v1 + v2 << std::endl;   
    std::cout << v1 + v2 << std::endl;
    std::cout << v1 - v2 << std::endl;
    std::cout << v1 * v2 << std::endl;
    std::cout << v1 * 2 << std::endl;
    std::cout << 2 * v1 << std::endl;
    
    return 0;
}