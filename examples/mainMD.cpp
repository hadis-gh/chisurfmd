#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <sstream>

#include <boost/program_options.hpp>

#include <adios2.h>

#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/Circle.h"
#include "chisurfmd/core/CircleDistribution.h"
#include "chisurfmd/core/ParticleDot.h"
#include "chisurfmd/core/ParticleOriented.h"
#include "chisurfmd/core/Utilities.h"
#include "chisurfmd/core/Energy.h"
#include "chisurfmd/core/FileIO.h"
#include "chisurfmd/md/Integration.h"
#include "chisurfmd/md/Thermostat.h"
#include "chisurfmd/potential/PotentialsU.h"

enum class ThermostatID {
    None = 0, VelocityScaling = 1, Berendsen = 2, NoseHoover = 3, Andersen = 4
};

#ifndef CHISURFMD_THERMOSTAT
#define CHISURFMD_THERMOSTAT ThermostatID::VelocityScaling
#endif

#ifndef CHISURFMD_PARTICLE
#define CHISURFMD_PARTICLE ParticleOriented
#endif

#ifndef CHISURFMD_POTENTIAL
#define CHISURFMD_POTENTIAL OrientedLJ
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace po = boost::program_options;
using Real = double;

using ParticleT = CHISURFMD_PARTICLE<Real>;
using Potential = CHISURFMD_POTENTIAL<ParticleT>;

void printOptions(const po::variables_map& vm) {
    std::cout << "\nCompile-time options:\n";
    std::cout << std::left << std::setw(20) << "Thermostat:" << TOSTRING(CHISURFMD_THERMOSTAT) << "\n"
              << std::setw(20) << "Particle type:" << TOSTRING(CHISURFMD_PARTICLE) << "\n"
              << std::setw(20) << "Potential type:" << TOSTRING(CHISURFMD_POTENTIAL) << "\n";

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
            } else if (value.type() == typeid(int)) {
                std::cout << option.second.as<int>();
            } else if (value.type() == typeid(std::vector<Real>)) {
                const auto& vec = option.second.as<std::vector<Real>>();
                std::cout << "[";
                for (int i = 0; i < (int)vec.size() - 1; i++) {
                    std::cout << vec[i] << ", ";
                }
                std::cout << vec.back() << "]";
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
        ("time,t",                po::value<Real>()->default_value(10.0),                     "max simulation time")
        ("dt",                    po::value<Real>()->default_value(.01),                      "integration step size")
        ("writeStateInterval",    po::value<Real>()->default_value(.05),                      "measurement State interval")
        ("writeEnergyInterval",   po::value<Real>()->default_value(.005),                     "measurement Energy interval")
        ("thermoInterval",        po::value<Real>()->default_value(.1),                       "interval after which to apply thermostat")
        ("temperature,T",         po::value<Real>()->default_value(.3),                       "temperature")
        ("particlesInit",         po::value<std::string>()->default_value("RANDOM"),          "particle initialization")
        ("chirality",             po::value<std::string>()->default_value("homochiral"),      "initial chirality: homochiral or racemic")
        ("alignment",             po::value<std::string>()->default_value("polar"),           "initial alignment: polar or apolar")        
        ("mass",                  po::value<Real>()->default_value(1.0),                      "mass of particles")       
        ("momentI",               po::value<Real>()->default_value(1.0),                      "moment of inersia")
        ("particleRadius",        po::value<Real>()->default_value(.5),                       "particle radius")
        ("seed",                  po::value<unsigned int>()->default_value(1),                "random seed")
        ("areaL",                 po::value<Real>()->default_value(20.0),                     "simulation size")
        ("fixRadius",             po::value<Real>()->default_value(20.0),                     "cut-off range for the dynamics neighbors")
        ("particlesDensity",      po::value<Real>(),                                          "packing density of particles")
        ("neighborDistances",     po::value<std::vector<Real>>()->multitoken()->default_value(std::vector<Real>{1.2, 2.0}, "1.2 2.0"),
                                                                                              "Distances for counting neighbors {x-y, omega}")
        ("particlesNum,n",        po::value<unsigned int>()->default_value(49),               "number of initial particles")
        ("saveFile",              po::value<std::string>()->default_value("outputs/run_0.bp"),"file path to save simulation output")
        ("relaxationTime",        po::value<Real>()->default_value(40.),                      "relaxation time for Berendsen thermostat")
        ("collisionFr",           po::value<Real>()->default_value(1.0 / 60.0),               "collision frequency for Andersen thermostat")
        ("integration",           po::value<std::string>()->default_value("VelocityVerlet"),  "integration method (velocity verlet/ euler)")
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

// ================================== simulation parameters ==================================
    
    const Real mass = vm["mass"].as<Real>();
    const Real momentI = vm["momentI"].as<Real>();
    const Real radius = vm["particleRadius"].as<Real>();
    unsigned int particlesNum = vm["particlesNum"].as<unsigned int>();
    Real areaL;

    if (vm.count("particlesDensity") > 0) {
        const Real density = vm["particlesDensity"].as<Real>();
        areaL = std::sqrt(particlesNum * M_PI * radius * radius / density);
    } else {
        areaL = vm["areaL"].as<Real>();
    }
    const Real boxPBC = vm["areaL"].as<Real>();
    const Real fixRadius = vm["fixRadius"].as<Real>();

    int speciesInd = 0;
    Species<Real> species1 {mass, momentI, radius};
    Species<Real> species2 {2.0 * mass, momentI, 0.5 * radius};
    std::vector<Species<Real>> allSpecies {species1, species2};

    std::mt19937 gen(vm["seed"].as<unsigned int>());

    auto particlesInit = vm["particlesInit"].as<std::string>(); 

    auto chirality = vm["chirality"].as<std::string>();
    auto alignment = vm["alignment"].as<std::string>();

    auto particles = initialParticles<ParticleT>(particlesNum, allSpecies, speciesInd, areaL, gen, particlesInit);
    
    if (particlesInit=="RANDOM") {
        assignParticleState(particles, chirality, alignment, gen);
    }
    
    particlesNum = particles.size();

    auto force = Potential::force(vm);
    auto potential = Potential::potential(vm);

    std::string method = vm["integration"].as<std::string>();
    auto integrationMethod = VelocityVerletStep<ParticleT, Potential::ForceType>;

    if      (method == "Euler")  integrationMethod = EulerStep<ParticleT, Potential::ForceType>;
    else if (method == "SEuler") integrationMethod = EulerSymplecticStep<ParticleT, Potential::ForceType>;
    else if (method != "VelocityVerlet") { std::cout << method << " integration wrong!"; return 1; }

    const std::vector<Real> neighborDistances = vm["neighborDistances"].as<std::vector<Real>>();
    const Real neighborCutoff = neighborDistances[0];

    const bool enableCapVelocity = vm["enableCapVelocity"].as<bool>();
    const std::vector<Real> maxVelocity = vm["maxVelocity"].as<std::vector<Real>>();

    const Real relaxationTime = vm["relaxationTime"].as<Real>();
    const Real collisionFrequency = vm["collisionFr"].as<Real>();
    const Real temperature = vm["temperature"].as<Real>();

    auto thermostat = [&]() {
        if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::None)
            return [](auto, auto) { return 1.; };
        else if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::VelocityScaling)
            return VelocityScalingThermostat<ParticleT>(temperature);
        else if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::Berendsen)
            return BerendsenThermostat<ParticleT>(vm["thermoInterval"].as<Real>(), temperature, relaxationTime);
        else if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::NoseHoover)
            return nullptr;
        else if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::Andersen)
            return AndersenThermostat<ParticleT>(vm["thermoInterval"].as<Real>(), collisionFrequency, temperature, gen);
    }();

// ================================== output files ==================================

    std::string adiosOutput = vm["saveFile"].as<std::string>();

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("SimulationOutput");

    constexpr int D = degreesOfFreedom<ParticleT>();

    adios2::Variable<Real> varT = io.DefineVariable<Real>("time");
    adios2::Variable<int8_t> varHandedness = io.DefineVariable<int8_t>("handedness", {particlesNum}, {0}, {particlesNum});
    adios2::Variable<int8_t> varAlignment = io.DefineVariable<int8_t>("alignment", {particlesNum}, {0}, {particlesNum});
    adios2::Variable<Real> varPositions = io.DefineVariable<Real>("positions", {particlesNum, D}, {0, 0}, {particlesNum, D});
    adios2::Variable<Real> varVelocities = io.DefineVariable<Real>("velocities", {particlesNum, D}, {0, 0}, {particlesNum, D});
    adios2::Variable<Real> varKineticEnergy = io.DefineVariable<Real>("kinetic energy", {1, 3}, {0, 0}, {1, 3});
    adios2::Variable<Real> varPotentialEnergy = io.DefineVariable<Real>("potential energy");
    adios2::Variable<Real> varNeighborCount = io.DefineVariable<Real>("number of neighbors", {1, neighborDistances.size()}, {0, 0}, {1, neighborDistances.size()});
    adios2::Variable<Real> varComVelocity = io.DefineVariable<Real>("center of mass velocity", {1, D}, {0, 0}, {1, D});
    adios2::Variable<Real> varComAngVelocity = io.DefineVariable<Real>("center of mass angular velocity");
    adios2::Variable<Real> varRealTemperature = io.DefineVariable<Real>("real temperature", {1, 3}, {0, 0}, {1, 3});
    adios2::Variable<Real> varOrientationOrder = io.DefineVariable<Real>("orientation order");

    io.DefineAttribute<std::string>("chirality", chirality);
    io.DefineAttribute<std::string>("alignment", alignment);
    io.DefineAttribute<Real>("particlesNum", particlesNum);
    io.DefineAttribute<Real>("temperature", temperature);
    io.DefineAttribute<Real>("radius", radius);
    io.DefineAttribute<Real>("mass", mass);
    io.DefineAttribute<Real>("momentI", momentI);
    io.DefineAttribute<Real>("areaL", areaL);
    io.DefineAttribute<Real>("neighborCutoff", neighborCutoff);
    
    adios2::Engine engine = io.Open(adiosOutput, adios2::Mode::Write);

// ================================== time intervals ==================================

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

    if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::None) {
        thermoStep = 2 * nsteps;
    }

    std::vector<int8_t> handednessVec(particlesNum);
    std::vector<int8_t> alignmentVec(particlesNum);
    std::vector<Real> positionsVec(particlesNum * D);
    std::vector<Real> velocitiesVec(particlesNum * D);

    clock_t startTime = clock();

// ================================== Integration Loop ==================================

    applyFixRadius(particles, fixRadius, areaL);
    
    while (step < nsteps) {
        engine.BeginStep();

        const auto nextEventStep = std::min({writeStateStep, writeEnergyStep, thermoStep, nsteps});
        Real currentTime = step * dt;
        integrate(particles, allSpecies, dt, (nextEventStep - step) * dt, boxPBC, force, integrationMethod);

        step = nextEventStep;

        engine.Put(varT, currentTime);

        if (enableCapVelocity) {
            capVelocity(particles, maxVelocity);
        }
        if (step == writeStateStep) {
            for (size_t j = 0; j < particles.size(); ++j) {
                const auto& p = particles[j];
                handednessVec[j] = p.handedness;
                alignmentVec[j]  = p.alignment;
                const auto pos = getGeneralizedPositions(p);
                const auto vel = getGeneralizedVelocities(p);
                for (int i = 0; i < D; ++i) {
                    positionsVec[j * D + i]  = pos[i];
                    velocitiesVec[j * D + i] = vel[i];
                }
            }

            engine.Put(varHandedness, handednessVec.data());
            engine.Put(varAlignment, alignmentVec.data());
            engine.Put(varPositions, positionsVec.data());
            engine.Put(varVelocities, velocitiesVec.data());
            writeStateStep = step + writeStateIntervalSteps;
        }
        if (step == writeEnergyStep) {
            auto kineticE = calKineticEnergy(particles, allSpecies);
            auto potentialE = calPotentialEnergy(particles, allSpecies, boxPBC, potential);

            auto neighborCount = calAveNeighborList(particles, neighborDistances, areaL);
            auto orientationOrder = calOrientationOrder(particles);

            auto COMvelocity = calCOMvelocity(particles, allSpecies);
            auto COMangularVelocity = calAngularMomentum2D(particles, allSpecies, boxPBC);
            
            auto realTemperature = calInternalTemperature(particles, allSpecies);

            engine.Put(varKineticEnergy, kineticE.data());
            engine.Put(varPotentialEnergy, potentialE);

            engine.Put(varNeighborCount, neighborCount.data());
            engine.Put(varOrientationOrder, orientationOrder);

            engine.Put(varComVelocity, COMvelocity.data());
            engine.Put(varComAngVelocity, COMangularVelocity);

            engine.Put(varRealTemperature, realTemperature.data());

            writeEnergyStep = step + writeEnergyIntervalSteps;
    
        }
        if constexpr (CHISURFMD_THERMOSTAT != ThermostatID::None) {
            if (step == thermoStep) {
                thermostat(particles, allSpecies);


                thermoStep = step + thermoIntervalSteps;
                if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::VelocityScaling) {
                    removeCOMVelocity(particles, allSpecies);
                }
            }
        }
        engine.EndStep();
    }
    engine.Close();
    
// ================================== run time output ==================================

    clock_t endTime = clock();

    Real timeTaken = Real(endTime - startTime) / CLOCKS_PER_SEC;

    int timeMin = static_cast<int>(timeTaken) / 60;
    double timeSec = std::fmod(timeTaken, 60);

    std::cout << "\nTime taken: " << timeMin << "m " << std::fixed << std::setprecision(2) << timeSec << "s\n";
    return 0;
}