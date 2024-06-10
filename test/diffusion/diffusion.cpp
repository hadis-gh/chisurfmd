#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include "lettuce/Vec.h"
#include "lettuce/Mat.h"

template <typename T>
struct Point {T x,y;};

std::vector <Vec<int, 2>> lattice_type(const std::string &name);
void simulate_randomwalk(std::vector<unsigned int>& count, std::vector<unsigned int>& endpoint_count, const int steps, const int samples);
void write_to_file(const std::vector<unsigned int>& count, const std::vector<unsigned int>& endpoint_count);
std::vector <Vec<int, 2>> directions = lattice_type("triangular");
const int L = 20;


int main() {
    std::vector<unsigned int> count(L*L, 0); 
    std::vector<unsigned int> endpoint_count(L*L, 0);

    int steps = 100;
    int traceStep = 10;
    int samples = 100000;

    simulate_randomwalk(count, endpoint_count, steps, samples);
    write_to_file(count, endpoint_count);

    return 0;
}


std::vector <Vec<int, 2>> lattice_type(const std::string &name){
    std::vector <Vec<int, 2>> directions;
    if (name=="triangular"){
        const Vec<int, 2> a1 ({1, 0});
        const Vec<int, 2> a2 ({0, 1});
        directions = {a1 , a2 , -a1 + a2, -a1, -a2, a1 -a2};
    }else if(name=="square"){
        const Vec<int, 2> a1 ({1, 0});
        const Vec<int, 2> a2 ({0, 1});
        directions = {a1 , a2 , -a1, -a2};        
    }
    return directions;
}

void simulate_randomwalk(std::vector<unsigned int>& count, std::vector<unsigned int>& endpoint_count, const int steps, const int samples){
    std::random_device rd;
    std::mt19937 gen(rd());     //previous one
    std::uniform_int_distribution<> distribution(0, directions.size()-1);

    for (int a=0; a<samples; ++a) {
    Vec<int, 2> current ({L/2, L/2});

    for (int i=0; i<steps; ++i) {
        int direction = distribution(gen);
        current += directions[direction];

        if (current >=0 && current< L)
            count.at(current[0] + current[1]*L) += 1;
    }
    if (current >=0 && current< L)
        endpoint_count.at(current[0] + current[1]*L) += 1;
    }
}

void write_to_file(const std::vector<unsigned int>& count, const std::vector<unsigned int>& endpoint_count) {
    std::ofstream output_file("atom_path.txt");
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int index = x + y * L;

            Vec<double, 2> lattice_coords({static_cast<double>(x), static_cast<double>(y)});
            Vec<double, 2> vec1({1, 0});
            Vec<double, 2> vec2({1 / 2, std::sqrt(3.0) / 2});
            Mat<double, 2, 2> basis({vec1, vec2});

            const auto cartesian = basis * lattice_coords;

            output_file << cartesian[0] << " " << cartesian[1] << " " 
                        << count[index] << " " << endpoint_count[index] << std::endl;
            // output_file << lattice_coords[0] << " " << lattice_coords[1] << " " 
            //             << count[index] << " " << endpoint_count[index] << std::endl;    
        }
    }
}

