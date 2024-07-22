#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <boost/program_options.hpp>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/Particle.h"
#include "lettuce/Utilities.h"
#include "lettuce/Integration.h"
#include "lettuce/Thermostat.h"
#include "lettuce/LennardJones.h"
#include "lettuce/ParticlesDistribution.h"
#include "lettuce/ParticlesIntersection.h"

enum class ThermostatID {
    None = 0, VelocityScaling = 1, Berendensen = 2, NoseHoover = 3, Andersen = 4
};

#ifndef LETTUCE_THERMOSTAT
#define LETTUCE_THERMOSTAT ThermostatID::None
#endif

namespace po = boost::program_options;
using Real = double;

int main(int argc, char* argv[]) {
    po::variables_map vm;
    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help,h", "print help")
        ("time,t",                po::value<double>()->default_value(100.0),        "max simulation time")
        ("dt",                    po::value<double>()->default_value(.001),         "integration step size")
        ("writeStateInterval",    po::value<unsigned int>()->default_value(1.),     "measurement State interval")
        ("writeEnergyInterval",   po::value<unsigned int>()->default_value(10.),     "measurement Energy interval")
        ("thermoInterval",        po::value<unsigned int>()->default_value(20.),     "interval after which to apply thermostat")
        ("temperature,T",         po::value<double>()->default_value(1.0),          "temperature")
        ("particleInit",          po::value<std::string>()->default_value("DLA"), "particle initialization")
        ("exclusionRadius",       po::value<double>()->default_value(.8),           "exclusion radius")
        ("seed",                  po::value<unsigned int>(),                        "random seed")
        ("areaL",                 po::value<double>()->default_value(100.0),        "simulation size")
        ("particlesNum,n",        po::value<int>()->default_value(10),              "number of initial particles");

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    if (vm.count("help") > 0) {
        std::cout << desc << std::endl;
        return 0;
    }

    int particlesNum = vm["particlesNum"].as<int>();

    auto gen = [&]() {
        if (vm.count("seed") > 0) {
            return std::mt19937(vm["seed"].as<unsigned int>());
        } else {
            std::random_device rd;
            return std::mt19937(rd());
        }
    }();

    const double sigma = 1;
    const double epsilon = 1;
    const double mass = 1;
    const double cutoff = 10;
    const double areaL = vm["areaL"].as<double>();
    const double radius = vm["exclusionRadius"].as<double>();

    const double dt = vm["dt"].as<double>();
    const double Time = vm["time"].as<double>();
    const unsigned int writeStateIntervalSteps = std::ceil(vm["writeStateInterval"].as<unsigned int>() / dt);
    const unsigned int writeEnergyIntervalSteps = std::ceil(vm["writeEnergyInterval"].as<unsigned int>() / dt);
    const unsigned int thermoIntervalSteps = std::ceil(vm["thermoInterval"].as<unsigned int>() / dt);

    const double relaxationTime = .1;
    const auto desiredTemperature = vm["temperature"].as<double>();

    LennardJonesForce<double> LJForce(epsilon, sigma, cutoff);
    LennardJonesPotential<double> LJPotential(epsilon, sigma, cutoff);

    auto Method = EulerStep<double, LennardJonesForce<double>>;

    Species<double> species1 {mass, radius};
    Species<double> species2 {2.0 * mass, 0.5 * radius};
    std::vector<Species<double>> allSpecies {species1, species2};
    int speciesInd = 1;

    std::vector<Particle<double>> particles = initialParticles(particlesNum, allSpecies, speciesInd, areaL, gen, vm["particleInit"].as<std::string>());

    std::ofstream positionFile("N_particle_PosMD.dat");
    std::ofstream kineticEnergyFile("N_particle_KineticEnergyMD.dat");
    std::ofstream PotentialEnergyFile("N_particle_PotentialEnergyMD.dat");

    auto thermostat = [&]() {
        if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
            return [](auto) { return 1.; };
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::VelocityScaling)
            return VelocityScalingThermostat<double>(desiredTemperature);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Berendensen)
            return BerendensenThermostat<double>(vm["thermoInterval"].as<unsigned int>() * dt, desiredTemperature, relaxationTime);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::NoseHoover)
            return nullptr;
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Andersen)
            return nullptr;
    }();

    const unsigned int nsteps = std::ceil(Time / dt);
    unsigned int step = 0;
    unsigned int writeStateStep = writeStateIntervalSteps;
    unsigned int writeEnergyStep = writeEnergyIntervalSteps;
    unsigned int thermoStep = thermoIntervalSteps;

    if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
        thermoStep = 2 * nsteps;

    clock_t startTime = clock();

    while (step < nsteps) {
        const auto nextEventStep = std::min({writeStateStep, writeEnergyStep, thermoStep, nsteps});
        integrate(particles, allSpecies, dt, (nextEventStep - step) * dt, LJForce, Method);
        step = nextEventStep;

        if (step == writeStateStep) {
            writePositionToFile(particles, positionFile);
            writeStateStep = step + writeStateIntervalSteps;
            std::cout << "written State at " << step * dt << " next " << writeStateStep * dt << std::endl;
        }
        if (step == writeEnergyStep) {
            writeKineticEToFile(particles, allSpecies, kineticEnergyFile);
            writePotentialEToFile(particles, allSpecies, LJPotential, PotentialEnergyFile);
            writeEnergyStep = step + writeEnergyIntervalSteps;
            std::cout << "written Energy at " << step * dt << " next " << writeEnergyStep * dt << std::endl;
        }        
        if constexpr (LETTUCE_THERMOSTAT != ThermostatID::None)
            if (step == thermoStep) {
                applyThermostat(particles, allSpecies, thermostat);
                thermoStep = step + thermoIntervalSteps;
                std::cout << "thermo at " << step * dt << " next " << thermoStep * dt << std::endl;
            }
    }
    clock_t endTime = clock();

    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    return 0;
}