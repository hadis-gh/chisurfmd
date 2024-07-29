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
#include "lettuce/CircleDistribution.h"


enum class ThermostatID {
    None = 0, VelocityScaling = 1, Berendsen = 2, NoseHoover = 3, Andersen = 4
};

#ifndef LETTUCE_THERMOSTAT
#define LETTUCE_THERMOSTAT ThermostatID::VelocityScaling
#endif

namespace po = boost::program_options;
using Real = double;

int main(int argc, char* argv[]) {
    po::variables_map vm;
    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help,h", "print help")
        ("time,t",                po::value<Real>()->default_value(20.0),            "max simulation time")
        ("dt",                    po::value<Real>()->default_value(.01),              "integration step size")
        ("writeStateInterval",    po::value<Real>()->default_value(1.),               "measurement State interval")
        ("writeEnergyInterval",   po::value<Real>()->default_value(1.),               "measurement Energy interval")
        ("thermoInterval",        po::value<Real>()->default_value(1.),              "interval after which to apply thermostat")
        ("temperature,T",         po::value<Real>()->default_value(.1),               "temperature")
        ("particleInit",          po::value<std::string>()->default_value("DLA"),     "particle initialization")
        ("exclusionRadius",       po::value<Real>()->default_value(.8),               "exclusion radius")
        ("seed",                  po::value<unsigned int>(),                          "random seed")
        ("areaL",                 po::value<Real>()->default_value(10.0),             "simulation size")
        ("particlesNum,n",        po::value<unsigned int>()->default_value(10),       "number of initial particles");

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    if (vm.count("help") > 0) {
        std::cout << desc << std::endl;
        return 0;
    }

    unsigned int particlesNum = vm["particlesNum"].as<unsigned int>();

    auto gen = [&]() {
        if (vm.count("seed") > 0) {
            return std::mt19937(vm["seed"].as<unsigned int>());
        } else {
            std::random_device rd;
            return std::mt19937(rd());
        }
    }();

    const Real sigma = 1;
    const Real epsilon = 1;
    const Real mass = 1;
    const Real cutoff = 10;
    const Real areaL = vm["areaL"].as<Real>();
    const Real radius = vm["exclusionRadius"].as<Real>();

    const Real dt = vm["dt"].as<Real>();
    const Real Time = vm["time"].as<Real>();

    const unsigned int writeStateIntervalSteps = std::ceil(vm["writeStateInterval"].as<Real>() / dt);
    const unsigned int writeEnergyIntervalSteps = std::ceil(vm["writeEnergyInterval"].as<Real>() / dt);
    const unsigned int thermoIntervalSteps = std::ceil(vm["thermoInterval"].as<Real>() / dt);

    const Real relaxationTime = 40.;
    const Real collisionFrequency = 1/40.;
    const auto desiredTemperature = vm["temperature"].as<Real>();

    LennardJonesForce<Real> LJForce(epsilon, sigma, cutoff);
    LennardJonesPotential<Real> LJPotential(epsilon, sigma, cutoff);

    auto Method = VelocityVerletStep<Real, LennardJonesForce<Real>>;

    Species<Real> species1 {mass, radius};
    Species<Real> species2 {2.0f * mass, 0.5f * radius};
    std::vector<Species<Real>> allSpecies {species1, species2};
    int speciesInd = 0;

    std::vector<Particle<Real>> particles = initialParticles(particlesNum, allSpecies, speciesInd, areaL, gen, vm["particleInit"].as<std::string>());
    
    std::ofstream initialParticles("initialConfigurationRead.dat");
    for (auto &p: particles){
        initialParticles << p.r[0] << " " << p.r[1] << " " << p.v[0] << " " << p.v[1] << std::endl;
    }

    std::ofstream positionFile("N_particle_PosMD.dat");
    std::ofstream kineticEnergyFile("N_particle_KineticEnergyMD.dat");
    std::ofstream PotentialEnergyFile("N_particle_PotentialEnergyMD.dat");

    auto thermostat = [&]() {
        if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
            return [](auto, auto) { return 1.; };
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::VelocityScaling)
            return VelocityScalingThermostat<Real>(desiredTemperature);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Berendsen)
            return BerendsenThermostat<Real>(vm["thermoInterval"].as<Real>(), desiredTemperature, relaxationTime);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::NoseHoover)
            return nullptr;
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Andersen){
            return AndersenThermostat<Real> (vm["thermoInterval"].as<Real>(), collisionFrequency, desiredTemperature, gen);
        }
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
        }
        if (step == writeEnergyStep) {
            writeKineticEToFile(particles, allSpecies, kineticEnergyFile);
            writePotentialEToFile(particles, allSpecies, LJPotential, PotentialEnergyFile);
            //write order parameter
            writeEnergyStep = step + writeEnergyIntervalSteps;
        }        
        if constexpr (LETTUCE_THERMOSTAT != ThermostatID::None)
            if (step == thermoStep) {
                thermostat(particles, allSpecies);
                thermoStep = step + thermoIntervalSteps;
            }    
    }

    std::ofstream configutation("/home/hadis/custom_vector/build/test/configuration.dat");
    for (auto p: particles){
        configutation << p.r[0] << " " << p.r[1] << " " << p.v[0] << " " << p.v[1] << std::endl;
    }

    std::ifstream configutation2("/home/hadis/custom_vector/build/test/configuration.dat");
    int numParticles = 0;
    std::string line;
    while (std::getline(configutation2, line)) {
        ++numParticles;
    }

    std::vector<Particle<Real>> particles2(numParticles);

    for (auto& p : particles2) {
        configutation2 >> p.r[0] >> p.r[1] >> p.v[0] >> p.v[1];
        std::cout << p.r[0] << " "<< p.r[1] << " " <<  p.v[0] << " "  << p.v[1] << std::endl;
    }


    clock_t endTime = clock();

    Real timeTaken = Real(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";    

    return 0;
}