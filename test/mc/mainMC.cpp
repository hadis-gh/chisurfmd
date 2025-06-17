#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <unordered_map>

#include <boost/program_options.hpp>
#include <adios2.h>

#include "lettuce/Vec.h"
#include "lettuce/Utilities.h"

namespace po = boost::program_options;
using Real = double;

void printOptions(const po::variables_map& vm) {
    if (vm["printOptions"].as<bool>()) {
        std::cout << "\nRuntime options:\n";
        for (const auto& option : vm) {
            std::cout << std::left << std::setw(20) << option.first << ": ";
            const auto& value = option.second.value();
            if (value.type() == typeid(std::string)) {
                std::cout << option.second.as<std::string>();
            } else if (value.type() == typeid(double)) {
                std::cout << option.second.as<double>();
            } else if (value.type() == typeid(unsigned int)) {
                std::cout << option.second.as<unsigned int>();
            } else if (value.type() == typeid(bool)) {
                std::cout << std::boolalpha << option.second.as<bool>();
            }
            std::cout << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help,h")
        ("steps",               po::value<unsigned int>()->default_value(100),                  "Steps per walker")
        ("particlesNum,n",      po::value<unsigned int>()->default_value(100),                  "Number of walkers")
        ("boxL",                po::value<Real>()->default_value(20.0),                         "Simulation box size (square)")
        ("latticeType",         po::value<std::string>()->default_value("square"),              "Lattice type")
        ("saveFile",            po::value<std::string>()->default_value("outputs/mc_0.bp"),     "ADIOS2 output file")
        ("seed",                po::value<unsigned int>(),                                      "Random seed")
        ("printOptions",        po::bool_switch()->default_value(true),                         "Print runtime options");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    printOptions(vm);

// ================================== parameters ==================================

    unsigned int steps = vm["steps"].as<unsigned int>();
    unsigned int N = vm["particlesNum"].as<unsigned int>();

    Real boxL = vm["boxL"].as<Real>();
    std::string latticeType = vm["latticeType"].as<std::string>();

    std::mt19937 gen = vm.count("seed") ? std::mt19937(vm["seed"].as<unsigned int>()) : std::mt19937(std::random_device{}());

    std::vector<Vec<int, 2>> directions = getLatticeDirections(latticeType);
    std::uniform_int_distribution<> dirDist(0, directions.size() - 1);

    int boxSize = static_cast<int>(boxL);
    std::vector<unsigned int> visitCounts(boxSize * boxSize, 0);
    std::vector<unsigned int> endpointCounts(boxSize * boxSize, 0);

// ================================== output ==================================

    std::string outputFile = vm["saveFile"].as<std::string>();
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("walk");

    auto varVisit = io.DefineVariable<unsigned int>("visitCounts", {boxSize * boxSize}, {0}, {boxSize * boxSize});
    auto varEnd = io.DefineVariable<unsigned int>("endpointCounts", {boxSize * boxSize}, {0}, {boxSize * boxSize});
    io.DefineAttribute<Real>("boxL", boxL);
    io.DefineAttribute<unsigned int>("steps", steps);
    io.DefineAttribute<unsigned int>("particlesNum", N);

    adios2::Engine engine = io.Open(outputFile, adios2::Mode::Write);

    // ================================== run simulation ==================================

    clock_t startTime = clock();

    for (unsigned int walker = 0; walker < N; ++walker) {
        Vec<int, 2> pos = {boxSize / 2, boxSize / 2};
        for (unsigned int s = 0; s < steps; ++s) {
            int dir = dirDist(gen);
            pos += directions[dir];

            if (pos[0] >= 0 && pos[0] < boxSize && pos[1] >= 0 && pos[1] < boxSize) {
                int idx = pos[0] + boxSize * pos[1];
                visitCounts[idx]++;
            }
        }
        if (pos[0] >= 0 && pos[0] < boxSize && pos[1] >= 0 && pos[1] < boxSize) {
            int idx = pos[0] + boxSize * pos[1];
            endpointCounts[idx]++;
        }
    }

    engine.BeginStep();
    engine.Put(varVisit, visitCounts.data());
    engine.Put(varEnd, endpointCounts.data());
    engine.EndStep();
    engine.Close();

// ================================== time ==================================
    clock_t endTime = clock();
    Real timeTaken = Real(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "\nTime taken: " << int(timeTaken / 60) << "m " << std::fixed << std::setprecision(2) << std::fmod(timeTaken, 60) << "s\n";

    return 0;
}
