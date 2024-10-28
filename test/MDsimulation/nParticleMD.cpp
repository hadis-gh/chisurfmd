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
#include "lettuce/ParticleOriented.h"
#include "lettuce/Utilities.h"
#include "lettuce/Integration.h"
#include "lettuce/Thermostat.h"
#include "lettuce/LennardJones.h"
#include "lettuce/LennardJonesOriented.h"
#include "lettuce/CircleDistribution.h"
#include "lettuce/LennarJonesParticles.h"

enum class ThermostatID {
    None = 0, VelocityScaling = 1, Berendsen = 2, NoseHoover = 3, Andersen = 4
};

#ifndef LETTUCE_THERMOSTAT
#define LETTUCE_THERMOSTAT ThermostatID::VelocityScaling
#endif

#ifndef LETTUCE_PARTICLE
#define LETTUCE_PARTICLE ParticleDot
#endif

#ifndef LETTUCE_POTENTIAL
#define LETTUCE_POTENTIAL LennardJones<ParticleT>
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace po = boost::program_options;
using Real = double;
// why not typedef double = Real; ?

// what is char* argv[] ? Does it something to do with lambda functions because of [] or just showing traditional lists?

using ParticleT = LETTUCE_PARTICLE<Real>;
using Potential = LETTUCE_POTENTIAL;

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
        ("seed",                  po::value<unsigned int>(),                          "random seed")
        ("areaL",                 po::value<Real>()->default_value(20.0),             "simulation size")
        ("particlesDensity",      po::value<Real>(),                                  "packing density of particles")
        ("neighborDist",          po::value<Real>()->default_value(1.2),              "Distance for counting neighbors")
        ("neighborDistances",     po::value<std::vector<Real>>(),                     "Distances for counting neighbors")
        ("particlesNum,n",        po::value<unsigned int>()->default_value(49),       "number of initial particles")
        ("saveParticles",         po::value<std::string>(),                           "file path to save final states")
        ("appendLog",             po::bool_switch(),                                  "append time series outputs")
        ("relaxationTime",        po::value<Real>()->default_value(40.),              "relaxation time for Berendsen thermostat")
        ("collisionFr",           po::value<Real>()->default_value(1.0 / 60.0),       "collision frequency for Andersen thermostat")
    ;

    Potential::initProgramOptions(desc);

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

    constexpr Real mass = 1;
    
    const Real radius = vm["exclusionRadius"].as<Real>();
    Real areaL = vm["areaL"].as<Real>();

    if (vm.count("particlesDensity") > 0) {
        const Real density = vm["particlesDensity"].as<Real>();
        const Real areaL = sqrt(particlesNum * M_PI * radius * radius / (density));
    }

    const Real boxPBC = vm["areaL"].as<Real>();
    const Real neighborDist = vm["neighborDist"].as<Real>();
    const std::vector<Real> neighborDistances {0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5};

    const Real relaxationTime = vm["relaxationTime"].as<Real>();
    const Real collisionFrequency = vm["collisionFr"].as<Real>();
    const Real desiredTemperature = vm["temperature"].as<Real>();

    auto force = Potential::force(vm);
    auto potential = Potential::potential(vm);
    // const auto cutoff = potential.cutoff();
    const auto& FORCE = Potential::force(vm);;

    auto integrationMethod = VelocityVerletStep<ParticleT, Potential::ForceType>;

    // auto integrationMethod = VelocityVerletStep<ParticleT, LennardJonesForce<Real>>;

    Species<Real> species1 {mass, radius};
    Species<Real> species2 {2.0f * mass, 0.5f * radius};
    std::vector<Species<Real>> allSpecies {species1, species2};
    int speciesInd = 0;
//encapsulation and abstraction of particle species instead of allSpecies + speciesInd + particles in the functions
    auto particlesInit = vm["particlesInit"].as<std::string>(); 

    auto particles = initialParticles<ParticleT>(particlesNum, allSpecies, speciesInd, areaL, gen, particlesInit);

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
    std::ofstream ComAngVelocityFile("N_particle_comAngVelocityMD.dat", openmode);
    std::ofstream TbeforeThermo("N_particle_TbeforeThermo.dat", openmode);
    std::ofstream TafterThermo("N_particle_TafterThermo.dat", openmode);
    std::ofstream realTbeforeThermo("N_particle_realTbeforeThermo.dat", openmode);
    std::ofstream realTafterThermo("N_particle_realTafterThermo.dat", openmode);

    std::cout << "\nthermostat: " << TOSTRING(LETTUCE_THERMOSTAT) << std::endl;
    std::cout << "parcticle type: " << TOSTRING(LETTUCE_PARTICLE) << std::endl;
    std::cout << "potential type: " << TOSTRING(LETTUCE_POTENTIAL) << std::endl;
    std::cout << "particle numbers: " << particlesNum << std::endl;
    std::cout << "particle initialization: " << particlesInit << "\n\n";

//print temperature
// same which prints force & potential / particle type (all run time and compile time parameter)
    auto thermostat = [&]() {
        if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
            return [](auto, auto) { return 1.; };
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::VelocityScaling)
            return VelocityScalingThermostat<ParticleT>(desiredTemperature);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Berendsen)
            return BerendsenThermostat<ParticleT>(vm["thermoInterval"].as<Real>(), desiredTemperature, relaxationTime);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::NoseHoover)
            return nullptr;
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Andersen) {
            return AndersenThermostat<ParticleT>(vm["thermoInterval"].as<Real>(), collisionFrequency, desiredTemperature, gen);
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
        integrate(particles, allSpecies, dt, (nextEventStep - step) * dt, boxPBC, force, integrationMethod);
        step = nextEventStep;

        if (step == writeStateStep) {
            writePositionToFile(particles, positionFile, dt, step);
            writeStateStep = step + writeStateIntervalSteps;
        }
        if (step == writeEnergyStep) {
            writeKineticEToFile(particles, allSpecies, kineticEnergyFile);
            writePotentialEToFile(particles, allSpecies, boxPBC, potential, PotentialEnergyFile);
            writeAverageNeighborToFile(particles, neighborDistances, NeighborCountFile, dt, step);
            writeComVelocityToFile(particles, allSpecies, ComVelocityFile, dt, step);
            writeComAngularVelocityToFile(particles, allSpecies, ComAngVelocityFile, dt, step);
            writeEnergyStep = step + writeEnergyIntervalSteps;
        }
        if constexpr (LETTUCE_THERMOSTAT != ThermostatID::None) {
            if (step == thermoStep) {
                writeTemperature(particles, allSpecies, TbeforeThermo);
                writeRealTemperature(particles, allSpecies, realTbeforeThermo);

                thermostat(particles, allSpecies);
                removeCOMVelocity(particles, allSpecies);
                removeCOMvelocityRotation2D(particles, allSpecies);
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
            auto q = getGeneralizedPositions(p);
            auto qDot = getGeneralizedVelocities(p);
            constexpr int D = degreesOfFreedom<ParticleT>();

            for (int a = 0; a < D; ++a) {
                configutation << std::setprecision(13) << std::scientific << q[a] << "\t";
            } for (int a = 0; a < D-1; ++a) {
                configutation << std::setprecision(13) << std::scientific << qDot[a] << "\t";
            }
            configutation << std::setprecision(13) << std::scientific << qDot[D-1] << std::endl;
        }
    }

    return 0;
}