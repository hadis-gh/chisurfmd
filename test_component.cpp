#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "Vec.h"
#include "Mat.h"

const double sqrt3 = std::sqrt(3.0);

template <typename T>
struct point{
    T x;
    T y;
};

int main()
{
    std::array <std::string, 2> lattice_types = {"square", "triangular"};
    std::string selected_lattice = "square";

    if(selected_lattice == "square"){
        Vec<int, 2> a1 ({1, 0});
        Vec<int, 2> a2 ({0, 1});

        std::array <Vec<int, 2>, 4> directions = {a1, a2, -a1, -a2};

    }else if(selected_lattice == "triangular"){
        Vec<double, 2> a1 ({1, 0});
        Vec<double, 2> a2 ({1/2, sqrt3/2});

        std::array <Vec<double, 2>, 6> directions = {a1 , a2 , -a1 + a2, -a1, -a2, a1 -a2};

    }else{
        std::cout << "Please choose the correct lattice!"<<std::endl;
    }

    int steps = 500;
    double probability = 0.5;
    double step_length = 1;

    point<double> origin ({0, 0});

    return 0;
}