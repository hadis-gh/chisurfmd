#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

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
#define LETTUCE_THERMOSTAT ThermostatID::None
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace po = boost::program_options;
using Real = double;
// why not typedef double = Real; ?

// what is char* argv[] ? Does it something to do with lambda functions because of [] or just showing traditional lists?

int main(int argc, char* argv[]) {
    po::variables_map vm;
    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help,h", "print help")
        ("time,t",                po::value<Real>()->default_value(10.0),             "max simulation time")
        ("dt",                    po::value<Real>()->default_value(.01),              "integration step size")
        ("writeStateInterval",    po::value<Real>()->default_value(.05),              "measurement State interval")
        ("writeEnergyInterval",   po::value<Real>()->default_value(.5),               "measurement Energy interval")
        ("thermoInterval",        po::value<Real>()->default_value(.1),               "interval after which to apply thermostat")
        ("temperature,T",         po::value<Real>()->default_value(.4),               "temperature")
        ("particlesInit",         po::value<std::string>()->default_value("RANDOM"),  "particle initialization")
        ("exclusionRadius",       po::value<Real>()->default_value(.8),               "exclusion radius")
        ("cutoff",                po::value<Real>()->default_value(10.),              "cutoff distance")
        ("seed",                  po::value<unsigned int>(),                          "random seed")
        ("areaL",                 po::value<Real>()->default_value(20.0),             "simulation size")
        ("particlesDensity",      po::value<Real>(),                                  "packing density of particles")
        ("neighborDist",          po::value<Real>()->default_value(1.2),              "Distance for counting neighbors")
        ("neighborDistances",     po::value<std::vector<Real>>(),                     "Distances for counting neighbors")
        ("particlesNum,n",        po::value<unsigned int>()->default_value(49),       "number of initial particles")
        ("saveParticles",         po::value<std::string>(),                           "file path to save final states")
        ("appendLog",             po::bool_switch(),                                  "append time series outputs")
        ("relaxationTime",        po::value<Real>()->default_value(40.),              "relaxation time for Berendsen thermostat")
        ("collisionFr",           po::value<Real>()->default_value(1.0 / 60.0),       "collision frequency for Andersen thermostat");

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

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

    constexpr Real sigma = 1;
    constexpr Real epsilon = 1;
    constexpr Real mass = 1;
// can't we ignore sigma and epsilon and consider them in units of r and u? (modified formula of LJ)
    
    const Real radius = vm["exclusionRadius"].as<Real>();
    Real areaL = vm["areaL"].as<Real>();

    if (vm.count("particlesDensity") > 0) {
        const Real density = vm["particlesDensity"].as<Real>();
        const Real areaL = sqrt(particlesNum * M_PI * radius * radius / (density));
    }

    const Real boxPBC = vm["areaL"].as<Real>();
    const Real cutoff = vm["cutoff"].as<Real>();
    const Real neighborDist = vm["neighborDist"].as<Real>();
    const std::vector<Real> neighborDistances {0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5};

    const Real relaxationTime = vm["relaxationTime"].as<Real>();
    const Real collisionFrequency = vm["collisionFr"].as<Real>();
    const Real desiredTemperature = vm["temperature"].as<Real>();

    LennardJonesForce<Real> LJForce(epsilon, sigma, cutoff);
    LennardJonesPotential<Real> LJPotential(epsilon, sigma, cutoff);

    auto integrationMethod = VelocityVerletStep<Real, LennardJonesForce<Real>>;

    Species<Real> species1 {mass, radius};
    Species<Real> species2 {2.0f * mass, 0.5f * radius};
    std::vector<Species<Real>> allSpecies {species1, species2};
    int speciesInd = 0;
//encapsulation and abstraction of particle species instead of allSpecies + speciesInd + particles in the functions

    std::vector<Particle<Real>> particles = initialParticles(particlesNum, allSpecies, speciesInd, areaL, gen, vm["particlesInit"].as<std::string>());
    writeInitialParticles(particles, radius);

    std::ios::openmode openmode = std::ios::trunc;
    if (vm["appendLog"].as<bool>()) {
        openmode = std::ios::app;
    }

    std::ofstream positionFile("N_particle_PosMD.dat", openmode);
    std::ofstream kineticEnergyFile("N_particle_KineticEnergyMD.dat", openmode);
    std::ofstream PotentialEnergyFile("N_particle_PotentialEnergyMD.dat", openmode);
    std::ofstream NeighborCountFile("N_particle_NeighborsMD.dat", openmode);
    std::ofstream ComVelocityFile("N_particle_comVelocityMD.dat", openmode);
    std::ofstream TbeforeThermo("N_particle_TbeforeThermo.dat", openmode);
    std::ofstream TafterThermo("N_particle_TafterThermo.dat", openmode);
    std::ofstream realTbeforeThermo("N_particle_realTbeforeThermo.dat", openmode);
    std::ofstream realTafterThermo("N_particle_realTafterThermo.dat", openmode);

    std::cout << "thermostat: " << TOSTRING(LETTUCE_THERMOSTAT) << std::endl;    

    auto thermostat = [&]() {
        if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
            return [](auto, auto) { return 1.; };
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::VelocityScaling)
            return VelocityScalingThermostat<Real>(desiredTemperature);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Berendsen)
            return BerendsenThermostat<Real>(vm["thermoInterval"].as<Real>(), desiredTemperature, relaxationTime);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::NoseHoover)
            return nullptr;
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Andersen) {
            return AndersenThermostat<Real>(vm["thermoInterval"].as<Real>(), collisionFrequency, desiredTemperature, gen);
        }
    }();

    const Real dt = vm["dt"].as<Real>();
    const Real Time = vm["time"].as<Real>();

    const size_t writeStateIntervalSteps = std::ceil(vm["writeStateInterval"].as<Real>() / dt);
    const size_t writeEnergyIntervalSteps = std::ceil(vm["writeEnergyInterval"].as<Real>() / dt);
    const size_t thermoIntervalSteps = std::ceil(vm["thermoInterval"].as<Real>() / dt);

    const size_t nsteps = std::ceil(Time / dt);

    size_t step = 0;
    size_t writeStateStep = writeStateIntervalSteps;
    size_t writeEnergyStep = writeEnergyIntervalSteps;
    size_t thermoStep = thermoIntervalSteps;

    if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None) {
        thermoStep = 2 * nsteps;
    }
    // considering line 184: if constexpr (LETTUCE_THERMOSTAT != ThermostatID::None) is not this redundent?

    clock_t startTime = clock();

    while (step < nsteps) {
        const auto nextEventStep = std::min({writeStateStep, writeEnergyStep, thermoStep, nsteps});
        integrate(particles, allSpecies, dt, (nextEventStep - step) * dt, boxPBC, LJForce, integrationMethod);
        step = nextEventStep;

        if (step == writeStateStep) {
            writePositionToFile(particles, positionFile, dt, step);
            writeStateStep = step + writeStateIntervalSteps;
        }
        if (step == writeEnergyStep) {
            writeKineticEToFile(particles, allSpecies, kineticEnergyFile);
            writePotentialEToFile(particles, allSpecies, boxPBC, LJPotential, PotentialEnergyFile);
            writeAverageNeighborToFile(particles, neighborDistances, NeighborCountFile, dt, step);
            writeComVelocityToFile(particles, allSpecies, ComVelocityFile, dt, step);
            writeEnergyStep = step + writeEnergyIntervalSteps;
        }
        if constexpr (LETTUCE_THERMOSTAT != ThermostatID::None) {
            if (step == thermoStep) {
                writeTemperature(particles, allSpecies, TbeforeThermo);
                writeRealTemperature(particles, allSpecies, realTbeforeThermo);

                thermostat(particles, allSpecies);
                thermoStep = step + thermoIntervalSteps;
                writeTemperature(particles, allSpecies, TafterThermo);
                writeRealTemperature(particles, allSpecies, realTafterThermo);
            }
        //deposition rate  (adding one particle to the list) -> make the option (by adding collisiotn frequency parameter) for running anderson thermostat after adding particle (but this one collistion frequency should be larger)
        }
    }

    clock_t endTime = clock();

    Real timeTaken = Real(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    if (vm.count("saveParticles") > 0) {
        std::ofstream configutation(vm["saveParticles"].as<std::string>());
        for (auto p: particles) {
            configutation << std::setprecision(13) << std::scientific << p.r[0] << "\t" << p.r[1] << "\t" << p.v[0] << "\t" << p.v[1] << std::endl;
        }
    }

    return 0;
}