#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

#include <boost/program_options.hpp>

#include <adios2.h>

#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"
#include "lettuce/Utilities.h"
#include "lettuce/Energy.h"
#include "lettuce/FileIO.h"
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
#define LETTUCE_PARTICLE ParticleOriented
#endif

#ifndef LETTUCE_POTENTIAL
#define LETTUCE_POTENTIAL LennardJones<ParticleT>
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace po = boost::program_options;
using Real = double;

using ParticleT = LETTUCE_PARTICLE<Real>;
using Potential = LETTUCE_POTENTIAL;

void printOptions(const po::variables_map& vm) {
    std::cout << "\nCompile-time options:\n";
    std::cout << std::left << std::setw(20) << "Thermostat:" << TOSTRING(LETTUCE_THERMOSTAT) << "\n"
              << std::setw(20) << "Particle type:" << TOSTRING(LETTUCE_PARTICLE) << "\n"
              << std::setw(20) << "Potential type:" << TOSTRING(LETTUCE_POTENTIAL) << "\n";

    if (vm["printOptions"].as<bool>()) {
        std::cout << "\nRuntime options:\n";
        for (const auto& option : vm) {
            const auto& value = option.second.value();
            std::cout << std::left << std::setw(20) << option.first << ": ";
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
        ("momentI",               po::value<Real>()->default_value(1.0),                      "moment of inersia")
        ("help,h", "print help")
        ("time,t",                po::value<Real>()->default_value(10.0),                     "max simulation time")
        ("dt",                    po::value<Real>()->default_value(.01),                      "integration step size")
        ("writeStateInterval",    po::value<Real>()->default_value(.05),                      "measurement State interval")
        ("writeEnergyInterval",   po::value<Real>()->default_value(.5),                       "measurement Energy interval")
        ("thermoInterval",        po::value<Real>()->default_value(.1),                       "interval after which to apply thermostat")
        ("temperature,T",         po::value<Real>()->default_value(.4),                       "temperature")
        ("particlesInit",         po::value<std::string>()->default_value("RANDOM"),          "particle initialization")
        ("exclusionRadius",       po::value<Real>()->default_value(.8),                       "exclusion radius")
        ("seed",                  po::value<unsigned int>(),                                  "random seed")
        ("areaL",                 po::value<Real>()->default_value(20.0),                     "simulation size")
        ("particlesDensity",      po::value<Real>(),                                          "packing density of particles")
        ("neighborDist",          po::value<Real>()->default_value(1.2),                      "Distance for counting neighbors")
        ("neighborDistances",     po::value<std::vector<Real>>(),                             "Distances for counting neighbors")
        ("particlesNum,n",        po::value<unsigned int>()->default_value(49),               "number of initial particles")
        ("saveParticles",         po::value<std::string>(),                                   "file path to save final states")
        ("appendLog",             po::bool_switch(),                                          "append time series outputs")
        ("relaxationTime",        po::value<Real>()->default_value(40.),                      "relaxation time for Berendsen thermostat")
        ("collisionFr",           po::value<Real>()->default_value(1.0 / 60.0),               "collision frequency for Andersen thermostat")
        ("integration",           po::value<std::string>()->default_value("VelocityVerlet"),  "integration method (velocity verlet/ euler)")
        ("saveAdios",         	  po::value<std::string>()->default_value("adios.bp"),        "file path to adios output file")
        ("enableCapVelocity",     po::bool_switch()->default_value(false),                    "Enable capping of velocities")
        ("maxVelocity",           po::value<std::vector<Real>>()->multitoken()->default_value(std::vector<Real>{1e5, 1e3}, "1e5 1e3"),
                                                                                              "capping amount for velocity {x-y, omega}")
        ("printOptions",          po::bool_switch()->default_value(true),                     "Print all runtime options")
;

    Potential::initProgramOptions(desc);

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help") > 0) {std::cout << desc << std::endl; return 0;}
    printOptions(vm);

    //simulation parameters
    constexpr Real mass = 1;
    const Real radius = vm["exclusionRadius"].as<Real>();
    unsigned int particlesNum = vm["particlesNum"].as<unsigned int>();
        if (vm.count("particlesDensity") > 0) {
        const Real density = vm["particlesDensity"].as<Real>();
        const Real areaL = sqrt(particlesNum * M_PI * radius * radius / (density));
    }
    Real areaL = vm["areaL"].as<Real>();
    const Real boxPBC = vm["areaL"].as<Real>();

    Species<Real> species1 {mass, vm["momentI"].as<Real>(), radius};
    Species<Real> species2 {2.0f * mass, vm["momentI"].as<Real>(), 0.5f * radius};
    std::vector<Species<Real>> allSpecies {species1, species2};
    int speciesInd = 0;
    auto particlesInit = vm["particlesInit"].as<std::string>(); 

    auto gen = [&]() {
        if (vm.count("seed") > 0) {
            return std::mt19937(vm["seed"].as<unsigned int>());
        } else {
            std::random_device rd;
            return std::mt19937(rd());
        }
    }();
    auto particles = initialParticles<ParticleT>(particlesNum, allSpecies, speciesInd, areaL, gen, particlesInit);
    writeInitialParticles(particles, radius);

    auto force = Potential::force(vm);
    auto potential = Potential::potential(vm);

    const Real neighborDist = vm["neighborDist"].as<Real>();
    const std::vector<Real> neighborDistances {0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5};
    
    const bool enableCapVelocity = vm["enableCapVelocity"].as<bool>();
    const std::vector<Real> maxVelocity = vm["maxVelocity"].as<std::vector<Real>>();

    const Real relaxationTime = vm["relaxationTime"].as<Real>();
    const Real collisionFrequency = vm["collisionFr"].as<Real>();
    const Real desiredTemperature = vm["temperature"].as<Real>();

    std::string method = vm["integration"].as<std::string>();
    auto integrationMethod = VelocityVerletStep<ParticleT, Potential::ForceType>;

    if (method=="VelocityVerlet") {integrationMethod = VelocityVerletStep<ParticleT, Potential::ForceType>;} 
    else if (method=="Euler")     {integrationMethod = EulerStep<ParticleT, Potential::ForceType>;} 
    else if (method=="SEuler")    {integrationMethod = EulerSymplecticStep<ParticleT, Potential::ForceType>;} 
    else {std::cout << method << " integration wrong!"; return 1;}

    auto thermostat = [&]() {
        if constexpr (LETTUCE_THERMOSTAT == ThermostatID::None)
            return [](auto, auto) { return 1.; };
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::VelocityScaling)
            return VelocityScalingThermostat<ParticleT>(desiredTemperature);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Berendsen)
            return BerendsenThermostat<ParticleT>(vm["thermoInterval"].as<Real>(), desiredTemperature, relaxationTime);
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::NoseHoover)
            return nullptr;
        else if constexpr (LETTUCE_THERMOSTAT == ThermostatID::Andersen)
            return AndersenThermostat<ParticleT>(vm["thermoInterval"].as<Real>(), collisionFrequency, desiredTemperature, gen);
    }();

    //output files
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("SimulationOutput");
    constexpr int D = degreesOfFreedom<ParticleT>();

    adios2::Variable<Real> varT = io.DefineVariable<Real>("time");
    adios2::Variable<Real> varPositions = io.DefineVariable<Real>("positions", {particlesNum, D}, {0, 0}, {particlesNum, D});
    adios2::Variable<Real> varVelocities = io.DefineVariable<Real>("velocities", {particlesNum, D}, {0, 0}, {particlesNum, D});
    adios2::Variable<Real> varKineticEnergy = io.DefineVariable<Real>("kinetic energy", {1, D}, {0, 0}, {1, D});
    adios2::Variable<Real> varPotentialEnergy = io.DefineVariable<Real>("potential energy");

    adios2::Engine engine = io.Open("SimulationOutput.bp", adios2::Mode::Write);

    std::ios::openmode openmode = std::ios::trunc;
    if (vm["appendLog"].as<bool>()) {openmode = std::ios::app;}

	// adios2::fstream oStream(vm["saveAdios"].as<std::string>(), adios2::fstream::out);

    // std::ofstream positionFile("N_particle_PosMD.dat", openmode);
    // std::ofstream kineticEnergyFile("N_particle_KineticEnergyMD.dat", openmode);
    // std::ofstream kEsAveFile("N_particle_KEsAveMD.dat", openmode);
    // std::ofstream PotentialEnergyFile("N_particle_PotentialEnergyMD.dat", openmode);
    // std::ofstream NeighborCountFile("N_particle_NeighborsMD.dat", openmode);
    // std::ofstream ComVelocityFile("N_particle_comVelocityMD.dat", openmode);
    // std::ofstream ComAngVelocityFile("N_particle_comAngVelocityMD.dat", openmode);
    // std::ofstream TbeforeThermo("N_particle_TbeforeThermo.dat", openmode);
    // std::ofstream TafterThermo("N_particle_TafterThermo.dat", openmode);
    // std::ofstream realTbeforeThermo("N_particle_realTbeforeThermo.dat", openmode);
    // std::ofstream realTafterThermo("N_particle_realTafterThermo.dat", openmode);

    //time intervals
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

    std::vector<Real> particlesVec(particlesNum * D);
    std::vector<Real> velocitiesVec(particlesNum * D);

    clock_t startTime = clock();
    engine.BeginStep();

    //main loop
    while (step < nsteps) {
        const auto nextEventStep = std::min({writeStateStep, writeEnergyStep, thermoStep, nsteps});
        integrate(particles, allSpecies, dt, (nextEventStep - step) * dt, boxPBC, force, integrationMethod);
        step = nextEventStep;
        engine.Put(varT, dt * step);
        if (enableCapVelocity) {
            capVelocity(particles, maxVelocity);
        }
        if (step == writeStateStep) {
            // writePositionToFile(particles, positionFile, dt, step);
            engine.Put(varPositions, particlesVec.data());
            engine.Put(varVelocities, particlesVec.data());
            writeStateStep = step + writeStateIntervalSteps;
        }
        if (step == writeEnergyStep) {
            // writeKineticEToFile(particles, allSpecies, oStream);
            // writeAverageKineticEnergies(particles, allSpecies, kEsAveFile);
            // writePotentialEToFile(particles, allSpecies, boxPBC, potential, PotentialEnergyFile);
            // writeAverageNeighborToFile(particles, neighborDistances, NeighborCountFile, dt, step);
            // writeComVelocityToFile(particles, allSpecies, ComVelocityFile, dt, step);
            // writeComAngularVelocityToFile(particles, allSpecies, ComAngVelocityFile, dt, step);
            auto kineticE = calInternalKineticEnergy(particles, allSpecies);
            auto potentialE = calPotentialEnergy(particles, allSpecies);
            engine.Put(varKineticEnergy, kineticE);
            engine.Put(varPotentialEnergy, potentialE);
            writeEnergyStep = step + writeEnergyIntervalSteps;
        }
        if constexpr (LETTUCE_THERMOSTAT != ThermostatID::None) {
            if (step == thermoStep) {
                // writeTemperature(particles, allSpecies, TbeforeThermo);
                // writeRealTemperature(particles, allSpecies, realTbeforeThermo);

                thermostat(particles, allSpecies);
                thermoStep = step + thermoIntervalSteps;
                if constexpr (LETTUCE_THERMOSTAT == ThermostatID::VelocityScaling) {
                    removeCOMVelocity(particles, allSpecies);
                    removeCOMvelocityRotation2D(particles, allSpecies, areaL);
                }
                
                // writeTemperature(particles, allSpecies, TafterThermo);
                // writeRealTemperature(particles, allSpecies, realTafterThermo);
            }
        //deposition rate  (adding one particle to the list) 
        //-> make the option (by adding collisiotn frequency parameter) for running anderson thermostat after adding particle (but this one collistion frequency should be larger)
        }

		// oStream.write<double>("time", dt * step, adios2::end_step);
		// oStream.write<double>("time", dt * step);
		// oStream.end_step();
    }
    engine.EndStep();
    engine.Close();
    
    clock_t endTime = clock();

    Real timeTaken = Real(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "\nTime taken: " << timeTaken << " seconds\n";

    if (vm.count("saveParticles") > 0) {
        std::ofstream configutation(vm["saveParticles"].as<std::string>());
        saveParticlesWithVelocities(configutation, particles);
    }

	// oStream.close();

    return 0;
}
