#include <iostream>
#include <fstream>
#include "lettuce/CircleDistribution.h"

int main() {
    int circles_number = 20;
    float length = 100.0;
    float radius = 1.0;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector <Circle<float>> circles = dist_circles (circles_number, length, radius, gen);
    for (auto c: circles)
        std::cout << c.c <<"\n";


    std::ofstream output_file("CirclesDistr.txt");
    for (auto c: circles)
        output_file << c.c <<"\n";

    return 0;
}