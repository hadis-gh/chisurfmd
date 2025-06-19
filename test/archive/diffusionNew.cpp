#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <sstream>

#include <boost/program_options.hpp>

#include <adios2.h>

#include "lettuce/core/Vec.h"
#include "lettuce/core/Circle.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"
#include "lettuce/core/Utilities.h"
#include "lettuce/core/Energy.h"
#include "lettuce/core/FileIO.h"
#include "lettuce/md/Integration.h"
#include "lettuce/md/Thermostat.h"
#include "lettuce/potential/IsotropicLJ.h"
#include "lettuce/potential/OrientedLJ.h"
#include "lettuce/potential/FieldCoupled.h"
#include "lettuce/potential/ChiralLJGeometric.h"
#include "lettuce/potential/PotentialFactory.h"
#include "lettuce/core/CircleDistribution.h"

#ifndef LETTUCE_PARTICLE
#define LETTUCE_PARTICLE ParticleDot
#endif

#ifndef LETTUCE_POTENTIAL
#define LETTUCE_POTENTIAL OrientedLJ
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace po = boost::program_options;
using Real = double;

using ParticleT = LETTUCE_PARTICLE<Real>;
using Potential = LETTUCE_POTENTIAL<ParticleT>;

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
    po::variables_map vm;
    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help,h", "print help")
        ("steps",               po::value<unsigned int>()->default_value(100),                  "Steps per walker")
        ("particlesNum,n",      po::value<unsigned int>()->default_value(100),                  "Number of walkers")
        ("boxL",                po::value<unsigned int>()->default_value(20),                   "Simulation box size (square)")
        ("latticeType",         po::value<std::string>()->default_value("square"),              "Lattice type")
        ("saveFile",            po::value<std::string>()->default_value("outputs/mc_0.bp"),     "ADIOS2 output file")
        ("seed",                po::value<unsigned int>(),                                      "Random seed")
        ("printOptions",        po::bool_switch()->default_value(true),                         "Print runtime options");

    Potential::initProgramOptions(desc);

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help") > 0) {std::cout << desc << std::endl; return 0;}
    printOptions(vm);

// ================================== parameters ==================================

    unsigned int steps = vm["steps"].as<unsigned int>();
    unsigned int N = vm["particlesNum"].as<unsigned int>();

    unsigned int boxL = vm["boxL"].as<unsigned int>();
    std::string latticeType = vm["latticeType"].as<std::string>();

    std::mt19937 gen = vm.count("seed") ? std::mt19937(vm["seed"].as<unsigned int>()) : std::mt19937(std::random_device{}());

    // std::vector<Vec<int, 2>> directions = getLatticeDirections(latticeType);
    std::vector<Vec<int, 2>> directions = {{{1, 0}}, {{-1, 0}}, {{0, 1}}, {{0, -1}}};

    std::uniform_int_distribution<> dirDist(0, directions.size() - 1);

    std::vector<unsigned int> visitCounts(boxL * boxL, 0);
    std::vector<unsigned int> endpointCounts(boxL * boxL, 0);

// ================================== output ==================================

    std::string outputFile = vm["saveFile"].as<std::string>();
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("walk");

    auto varVisit = io.DefineVariable<unsigned int>("visitCounts", {boxL * boxL}, {0}, {boxL * boxL});
    auto varEnd = io.DefineVariable<unsigned int>("endpointCounts", {boxL * boxL}, {0}, {boxL * boxL});
    io.DefineAttribute<unsigned int>("boxL", boxL);
    io.DefineAttribute<unsigned int>("steps", steps);
    io.DefineAttribute<unsigned int>("particlesNum", N);

    adios2::Engine engine = io.Open(outputFile, adios2::Mode::Write);

    // ================================== run simulation ==================================

    clock_t startTime = clock();

    for (unsigned int walker = 0; walker < N; ++walker) {
        Vec<int, 2> pos = {{static_cast<int>(boxL / 2), static_cast<int>(boxL / 2)}};
        
        for (unsigned int s = 0; s < steps; ++s) {
            int dir = dirDist(gen);
            pos += directions[dir];

            if (pos[0] >= 0 && pos[0] < boxL && pos[1] >= 0 && pos[1] < boxL) {
                int idx = pos[0] + boxL * pos[1];
                visitCounts[idx]++;
            }
        }
        if (pos[0] >= 0 && pos[0] < boxL && pos[1] >= 0 && pos[1] < boxL) {
            int idx = pos[0] + boxL * pos[1];
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
