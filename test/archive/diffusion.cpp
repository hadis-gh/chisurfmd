#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <chrono>

#include "lettuce/core/Vec.h"
#include "lettuce/core/Mat.h"

void simulateRandomWalk(std::vector<unsigned int>& visitCounts,
                        std::vector<unsigned int>& endpointCounts,
                        const int numSteps, const int numWalkers,
                        const int boxSize, const std::string& latticeType);

void writeToFile(const std::vector<unsigned int>& visitCounts,
                 const std::vector<unsigned int>& endpointCounts,
                 const int boxSize);

std::vector<Vec<int, 2>> getLatticeDirections(const std::string& latticeType);

int main() {
    const int boxSize = 20;
    std::vector<unsigned int> visitCountsGrid(boxSize * boxSize, 0);
    std::vector<unsigned int> endpointCountsGrid(boxSize * boxSize, 0);

    const int numSteps = 100;
    const int numWalkers = 100000;
    const std::string latticeType = "triangular";  // or "square"

    simulateRandomWalk(visitCountsGrid, endpointCountsGrid, numSteps, numWalkers, boxSize, latticeType);
    writeToFile(visitCountsGrid, endpointCountsGrid, boxSize);

    return 0;
}

void simulateRandomWalk(std::vector<unsigned int>& visitCounts,
                        std::vector<unsigned int>& endpointCounts,
                        const int numSteps, const int numWalkers,
                        const int boxSize, const std::string& latticeType) {

    std::vector<Vec<int, 2>> movementDirections = getLatticeDirections(latticeType);

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<> directionDist(0, movementDirections.size() - 1);

    for (int walker = 0; walker < numWalkers; ++walker) {
        Vec<int, 2> position({boxSize / 2, boxSize / 2});

        for (int step = 0; step < numSteps; ++step) {
            int dir = directionDist(rng);
            position += movementDirections[dir];

            if (position[0] >= 0 && position[0] < boxSize &&
                position[1] >= 0 && position[1] < boxSize) {
                int index = position[0] + position[1] * boxSize;
                visitCounts[index]++;
            }
        }

        if (position[0] >= 0 && position[0] < boxSize &&
            position[1] >= 0 && position[1] < boxSize) {
            int index = position[0] + position[1] * boxSize;
            endpointCounts[index]++;
        }
    }
}

void writeToFile(const std::vector<unsigned int>& visitCounts,
                 const std::vector<unsigned int>& endpointCounts,
                 const int boxSize) {
    std::ofstream output("atomPath.txt");

    for (int y = 0; y < boxSize; ++y) {
        for (int x = 0; x < boxSize; ++x) {
            int index = x + y * boxSize;

            Vec<double, 2> latticeCoords({static_cast<double>(x), static_cast<double>(y)});
            Vec<double, 2> a1({1.0, 0.0});
            Vec<double, 2> a2({0.5, std::sqrt(3.0) / 2.0});
            Mat<double, 2, 2> basis({a1, a2});
            Vec<double, 2> cartesianCoords = basis * latticeCoords;

            output << cartesianCoords[0] << " " << cartesianCoords[1] << " "
                   << visitCounts[index] << " " << endpointCounts[index] << "\n";
        }
    }
}

std::vector<Vec<int, 2>> getLatticeDirections(const std::string& latticeType) {
    if (latticeType == "triangular") {
        Vec<int, 2> a1({1, 0});
        Vec<int, 2> a2({0, 1});
        return {a1, a2, -a1 + a2, -a1, -a2, a1 - a2};
    } else if (latticeType == "square") {
        Vec<int, 2> a1({1, 0});
        Vec<int, 2> a2({0, 1});
        return {a1, a2, -a1, -a2};
    } else {
        throw std::invalid_argument( "invalid lattice type" );
    }
}
