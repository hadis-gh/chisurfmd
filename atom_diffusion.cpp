#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <random>
#include "Vec.h"
#include "Mat.h"

template <class T>
struct Point {
    T x;
    T y;
};

const double sqrt3 = std::sqrt(3.0);

const Vec<int, 2> a1 ({1, 0});
const Vec<int, 2> a2 ({0, 1});

const std::array <Vec<int, 2>, 6> directions = {a1 , a2 , -a1 + a2, -a1, -a2, a1 -a2};

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, directions.size()-1);

    const int min = 0;
    const int L = 20;
    std::vector<unsigned int> hist(L*L, 0); 
    std::vector<unsigned int> hist100(L*L, 0);

    int steps = 100;
    int traceStep = 10;
    int samples = 10000000;

    for (int a=0; a<samples; ++a) {
        Vec<int, 2> current ({L/2, L/2});

        for (int i=0; i<steps; ++i) {
            int direction = dis(gen);
            current += directions[direction];

            if (current >=0 && current< L)
                hist.at(current[0] + current[1]*L) += 1;
        }
        if (current >=0 && current< L)
            hist100.at(current[0] + current[1]*L) += 1;
    }

    Vec<double, 2> vec1({1, 0});
    Vec<double, 2> vec2({1/2, sqrt3/2});
    Mat<double, 2, 2> basis({vec1, vec2});

    std::ofstream outFile("atom_path.txt");
    for (int y = min; y< min + L; ++y){
        for (int x = min; x<min + L; ++x){
            const int i = x-min + (y-min)*L;
            Vec<double, 2> red_path ({double(x),double(y)});
            const auto cart = basis * red_path;
            outFile << cart[0] << " " << cart[1] << " " << hist[i]<< " "<< hist100[i]<< std::endl;
        }
    }

    return 0;
}