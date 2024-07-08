#include <iostream>
#include <fstream>
#include "lettuce/CircleDistribution.h"

int main() {
    int circlesNumber = 40;
    float length = 100.0;
    float radius = 2.0;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<Circle<float>> circles = distRandomCircles(circlesNumber, length, radius, gen);

    std::ofstream output_file("CirclesDistr.txt");
    for (auto c: circles)
        output_file << c.c <<"\n";

    return 0;
}