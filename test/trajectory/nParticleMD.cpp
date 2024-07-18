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
#include "lettuce/CommandLineOptions.h"

enum class ThermostatID
{
    None = 0, Berendensen = 1, Andersen, VelocityScaling, NoseHoover
};
#ifndef LETTUCE_THERMOSTAT
#define LETTUCE_THERMOSTAT ThermostatID::None
#endif

namespace po = boost::program_options;
using Real = double;

int main(int argc, char* argv[]){

    po::variables_map vm = parseCommandLineOptions(argc, argv);

    int particlesNum = vm["particlesNum"].as<int>();

    auto gen = [&]() {
        if(vm.count("seed") > 0)
        {
            return std::mt19937(vm["seed"].as<unsigned int>());
        }
        else
        {
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
    const unsigned int writeInterval = std::ceil(vm["writeInterval"].as<double>()/dt);
    const unsigned int thermoInterval = std::ceil(vm["thermoInterval"].as<double>()/dt);

    const double relaxationTime = .1;
    const auto desiredTemperature = vm["temperature"].as<double>();

    LennardJonesForce<double> LJForce(epsilon, sigma, cutoff);
    LennardJonesPotential<double> LJPotential(epsilon, sigma, cutoff);

    auto Method = EulerStep<double, LennardJonesForce<double>>;

    Species<double> species1 {mass, radius};
    Species<double> species2 {2.0 * mass, 0.5 * radius};
    std::vector <Species<double>> allSpecies {species1, species2};
    int speciesInd = 1;

    std::vector<Particle<double>> particles = initialParticles (particlesNum, allSpecies, speciesInd, areaL, gen, vm["particleInit"].as<std::string>());

    std::ofstream positionFile ("N_particle_PosMD.dat");
    std::ofstream kineticEnergyFile ("N_particle_KineticEnergyMD.dat");
    std::ofstream PotentialEnergyFile ("N_particle_PotentialEnergyMD.dat");


    auto thermostat = [&](){
        if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
            return [](auto){ return 1.; };
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Berendensen)
            return Berendensen<double>(vm["thermoInterval"].as<double>(), desiredTemperature, relaxationTime);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Andersen)
            return 0; //AndersenThermostat<double>(desiredTemperature);
    }();

    unsigned int time = 0;
    unsigned int writeTime = time, thermoTime = thermoInterval;
    const unsigned int nsteps = std::ceil(Time / dt);

    if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
        thermoTime = 2 * nsteps;


    clock_t startTime = clock();

    while (time < nsteps){
        const auto nextEventTime = std::min(std::min(writeTime, thermoTime), nsteps);
        integrate(particles, allSpecies, dt, (nextEventTime - time)*dt, LJForce, Method);
        time = nextEventTime;

        if(time == writeTime){
            writePositionToFile(particles, positionFile);
            writeKineticEToFile(particles, allSpecies, kineticEnergyFile);
            writePotentialEToFile(particles, allSpecies, LJPotential, PotentialEnergyFile);
            writeTime = time + writeInterval;
            std::cout << "written at " << time*dt << " next " << writeTime*dt << std::endl;
        }
        if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
            if(time == thermoTime){
                applyThermostat(particles, allSpecies, thermostat);
                thermoTime = time + thermoInterval;
                std::cout << "thermo at " << time*dt << " next " << thermoTime*dt << std::endl;
            }
    }

    clock_t endTime = clock();

    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    return 0;
}