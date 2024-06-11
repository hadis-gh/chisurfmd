#include <iostream>
#include <fstream>
#include "lettuce/CircleDistribution.h"

int main() {
    int circles_number = 40;
    float length = 100.0;
    float radius = 2.0;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::ofstream output_file("CirclesDistr.txt");
    for (auto c: circles)
        output_file << c.c <<"\n";

    return 0;
}