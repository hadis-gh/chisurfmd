#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <sstream>
#include <typeinfo>

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
#define CHISURFMD_PARTICLE ParticleDot
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
            }
            std::cout << std::endl;
        }
    }
}

template<typename P>
void assignPhiIfOriented(P& particle, std::mt19937& gen) {
    if constexpr (degreesOfFreedom<P>() > 2) {
        std::uniform_real_distribution<Real> phiDist(0.0, 2.0 * M_PI);
        particle.phi = phiDist(gen);
    }
}

int main(int argc, char* argv[]) {
    po::variables_map vm;
    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help,h", "print help")
        ("time,t",                po::value<Real>()->default_value(30.0),                     "max simulation time")
        ("dt",                    po::value<Real>()->default_value(.00005),                   "integration step size")
        ("writeStateInterval",    po::value<Real>()->default_value(.01),                      "measurement State interval")
        ("writeEnergyInterval",   po::value<Real>()->default_value(.1),                       "measurement Energy interval")
        ("thermoInterval",        po::value<Real>()->default_value(.05),                      "interval after which to apply thermostat")
        ("temperature,T",         po::value<Real>()->default_value(.3),                       "temperature")
        ("particlesInit",         po::value<std::string>()->default_value("TWO"),             "particle initialization")
        ("chirality",             po::value<std::string>()->default_value("homochiral"),      "initial chirality: homochiral or racemic")
        ("alignment",             po::value<std::string>()->default_value("polar"),           "initial alignment: polar or apolar")   
        ("mass",                  po::value<Real>()->default_value(1.0),                      "mass of particles")       
        ("momentI",               po::value<Real>()->default_value(1.0),                      "moment of inersia")
        ("particleRadius",       po::value<Real>()->default_value(.8),                        "particle radius")
        ("seed",                  po::value<unsigned int>(),                                  "random seed")
        ("areaL",                 po::value<Real>()->default_value(50.0),                     "simulation size")
        ("particlesDensity",      po::value<Real>(),                                          "packing density of particles")
        ("neighborDistances",     po::value<std::vector<Real>>()->multitoken()->default_value(std::vector<Real>{1.2, 2.0}, "1.2 2.0"),
                                                                                              "Distances for counting neighbors {x-y, omega}")
        ("particlesNum,n",        po::value<unsigned int>()->default_value(2),                "number of initial particles")
        ("particlesNumMax",       po::value<unsigned int>()->default_value(49),               "number of final particles")
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
    unsigned int particlesNumMax = vm["particlesNumMax"].as<unsigned int>();
        if (vm.count("particlesDensity") > 0) {
        const Real density = vm["particlesDensity"].as<Real>();
        const Real areaL = sqrt(particlesNum * M_PI * radius * radius / (density));
    }
    Real areaL = vm["areaL"].as<Real>();
    const Real boxPBC = vm["areaL"].as<Real>();

    int speciesInd = 0;
    Species<Real> species1 {mass, momentI, radius};
    Species<Real> species2 {2.0f * mass, momentI, 0.5f * radius};
    std::vector<Species<Real>> allSpecies {species1, species2};

    auto particlesInit = vm["particlesInit"].as<std::string>(); 
    
    auto chirality = vm["chirality"].as<std::string>();
    auto alignment = vm["alignment"].as<std::string>();

    auto gen = [&]() {
        if (vm.count("seed") > 0) {
            return std::mt19937(vm["seed"].as<unsigned int>());
        } else {
            std::random_device rd;
            return std::mt19937(rd());
        }
    }();

    auto force = Potential::force(vm);
    auto potential = Potential::potential(vm);

    std::string method = vm["integration"].as<std::string>();
    auto integrationMethod = VelocityVerletStep<ParticleT, Potential::ForceType>;

    if (method=="VelocityVerlet") {integrationMethod = VelocityVerletStep<ParticleT, Potential::ForceType>;} 
    else if (method=="Euler")     {integrationMethod = EulerStep<ParticleT, Potential::ForceType>;} 
    else if (method=="SEuler")    {integrationMethod = EulerSymplecticStep<ParticleT, Potential::ForceType>;} 
    else {std::cout << method << " integration wrong!"; return 1;}

    const std::vector<Real> neighborDistances = vm["neighborDistances"].as<std::vector<Real>>();
    const Real neighborCutoff = neighborDistances[0];

    const bool enableCapVelocity = vm["enableCapVelocity"].as<bool>();
    const std::vector<Real> maxVelocity = vm["maxVelocity"].as<std::vector<Real>>();

    const Real relaxationTime = vm["relaxationTime"].as<Real>();
    const Real collisionFrequency = vm["collisionFr"].as<Real>();

    const Real collisionFrDeposit = vm["collisionFr"].as<Real>();
    
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
    adios2::Variable<Real> varPositions = io.DefineVariable<Real>("positions", {particlesNumMax, D}, {0, 0}, {particlesNumMax, D});
    adios2::Variable<Real> varVelocities = io.DefineVariable<Real>("velocities", {particlesNumMax, D}, {0, 0}, {particlesNumMax, D});
    adios2::Variable<Real> varKineticEnergy = io.DefineVariable<Real>("kinetic energy", {1, 3}, {0, 0}, {1, 3});
    adios2::Variable<Real> varPotentialEnergy = io.DefineVariable<Real>("potential energy");
    adios2::Variable<Real> varNeighborCount = io.DefineVariable<Real>("number of neighbors", {1, neighborDistances.size()}, {0, 0}, {1, neighborDistances.size()});
    adios2::Variable<Real> varComVelocity = io.DefineVariable<Real>("center of mass velocity", {1, D}, {0, 0}, {1, D});
    adios2::Variable<Real> varComAngVelocity = io.DefineVariable<Real>("center of mass angular velocity");
    adios2::Variable<Real> varRealTemperature = io.DefineVariable<Real>("real temperature", {1, 3}, {0, 0}, {1, 3});
    adios2::Variable<Real> varOrientationOrder = io.DefineVariable<Real>("orientation order");

    io.DefineAttribute<std::string>("chirality", chirality);
    io.DefineAttribute<std::string>("alignment", alignment);
    io.DefineAttribute<Real>("temperature", temperature);
    io.DefineAttribute<Real>("radius", radius);
    io.DefineAttribute<Real>("mass", mass);
    io.DefineAttribute<Real>("momentI", momentI);
    io.DefineAttribute<Real>("areaL", areaL);
    io.DefineAttribute<Real>("particlesNum", particlesNum);
    io.DefineAttribute<Real>("particlesNumMax", particlesNumMax);
    io.DefineAttribute<Real>("neighborCutoff", neighborCutoff);

    const int8_t INVALID_VALUE = std::numeric_limits<int8_t>::min();

    std::vector<int8_t> handednessVec(particlesNumMax, INVALID_VALUE);
    std::vector<int8_t> alignmentVec(particlesNumMax, INVALID_VALUE);
    std::vector<Real> positionsVec(particlesNumMax * D, INVALID_VALUE);
    std::vector<Real> velocitiesVec(particlesNumMax * D, INVALID_VALUE);
    
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
        
    clock_t startTime = clock();

// ================================== Deposition Loop ==================================
    
    auto particles = initialParticles<ParticleT>(particlesNum, allSpecies, speciesInd, areaL, gen, particlesInit);
    
    auto prevSize = particles.size();
    std::cout << "\n _______ system initial size: " << particles.size() << " _______ \n" << std::endl;

    moveParticlesToCenter(particles, allSpecies, areaL);

    for (int a=prevSize;  particles.size()<particlesNumMax; ++a) {
        auto newPos = depositeDLAinfo3(particles, allSpecies, speciesInd, areaL, gen);
        
        ParticleT newParticle(speciesInd, newPos);

        assignPhiIfOriented(newParticle, gen);
        assignNewParticleState(newParticle, particles, chirality, alignment, gen);

        particles.push_back(newParticle);
        std::cout << "\rNew particle added! System size: " << particles.size() << std::flush;

        resetVelocitiesRandom(particles, allSpecies, temperature, gen);

        step = 0;
        writeStateStep = writeStateIntervalSteps;
        writeEnergyStep = writeEnergyIntervalSteps;
        thermoStep = thermoIntervalSteps;

// ================================== Integration Loop ==================================
        
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
                handednessVec.resize(particlesNumMax, INVALID_VALUE);
                alignmentVec.resize(particlesNumMax, INVALID_VALUE);
                positionsVec.resize(particlesNumMax * D, INVALID_VALUE);
                velocitiesVec.resize(particlesNumMax * D, INVALID_VALUE);
    
                for (size_t i = 0; i < particles.size(); ++i) {
                    auto pos = getGeneralizedPositions(particles[i]);
                    auto vel = getGeneralizedVelocities(particles[i]);
                    for (int j = 0; j < D; ++j) {
                        positionsVec[i * D + j] = pos[j];
                        velocitiesVec[i * D + j] = vel[j];
                    }
                    handednessVec[i] = particles[i].handedness;
                    alignmentVec[i] = particles[i].alignment;
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
    
                engine.Put(varKineticEnergy, kineticE.data());
                engine.Put(varPotentialEnergy, potentialE);
    
                engine.Put(varNeighborCount, neighborCount.data());
                engine.Put(varOrientationOrder, orientationOrder);

                engine.Put(varComVelocity, COMvelocity.data());
                engine.Put(varComAngVelocity, COMangularVelocity);
    
                writeEnergyStep = step + writeEnergyIntervalSteps;
            }
            if constexpr (CHISURFMD_THERMOSTAT != ThermostatID::None) {
                if (step == thermoStep) {
                    thermostat(particles, allSpecies);
    
                    auto realTemperature = calInternalTemperature(particles, allSpecies);
    
                    engine.Put(varRealTemperature, realTemperature.data());
    
                    thermoStep = step + thermoIntervalSteps;
                    if constexpr (CHISURFMD_THERMOSTAT == ThermostatID::VelocityScaling) {
                        removeCOMVelocity(particles, allSpecies);
                    }
                }
            }
            engine.EndStep();
        }    
    }

    engine.Close();

    // ================================== run time output ==================================

    clock_t endTime = clock();

    Real timeTaken = Real(endTime - startTime) / CLOCKS_PER_SEC;

    int timeMin = static_cast<int>(timeTaken) / 60;
    double timeSec = std::fmod(timeTaken, 60);
    
    std::cout << "\nDone!" << std::endl;
    std::cout << "\nTime taken: " << timeMin << "m " << std::fixed << std::setprecision(2) << timeSec << "s\n";

    return 0;
}
