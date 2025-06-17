#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <sstream>

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
#include "lettuce/IsotropicLJ.h"
#include "lettuce/OrientedLJ.h"
#include "lettuce/FieldCoupled.h"
#include "lettuce/ChiralLJGeometric.h"
#include "lettuce/PotentialFactory.h"
#include "lettuce/CircleDistribution.h"


#ifndef LETTUCE_PARTICLE
#define LETTUCE_PARTICLE ParticleOriented
#endif

#ifndef LETTUCE_POTENTIAL
#define LETTUCE_POTENTIAL ChiralLJGeometric<ParticleT>
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace po = boost::program_options;
using Real = double;

using ParticleT = LETTUCE_PARTICLE<Real>;
using Potential = LETTUCE_POTENTIAL;

void printOptions(const po::variables_map& vm) {
    std::cout << "\nCompile-time options:\n";
    std::cout << std::setw(20) << "Particle type:" << TOSTRING(LETTUCE_PARTICLE) << "\n"
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
        ("help,h", "print help")
        ("steps",                 po::value<unsigned int>()->default_value(100),              "total steps of simulation")
        ("particlesNum,n",        po::value<unsigned int>()->default_value(49),               "number of initial particles")
        ("particlesInit",         po::value<std::string>()->default_value("RANDOM"),          "particle initialization")
        ("radius",                po::value<Real>()->default_value(.5),                       "partilce's radius")
        ("mass",                  po::value<Real>()->default_value(1.0),                      "mass of particles")       
        ("momentI",               po::value<Real>()->default_value(1.0),                      "moment of inersia")
        ("seed",                  po::value<unsigned int>(),                                  "random seed")
        ("areaL",                 po::value<Real>()->default_value(20.0),                     "simulation size")
        ("saveFile",              po::value<std::string>()->default_value("outputs/run_0.bp"),"file path to save simulation output")
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
    const Real radius = vm["radius"].as<Real>();
    unsigned int particlesNum = vm["particlesNum"].as<unsigned int>();

    Real areaL = vm["areaL"].as<Real>();
    const Real boxPBC = vm["areaL"].as<Real>();

    int speciesInd = 0;
    Species<Real> species1 {mass, momentI, radius};
    Species<Real> species2 {2.0f * mass, momentI, 0.5f * radius};
    std::vector<Species<Real>> allSpecies {species1, species2};

    auto particlesInit = vm["particlesInit"].as<std::string>(); 
    auto particlesType = vm["particlesType"].as<std::string>(); 

    auto gen = [&]() {
        if (vm.count("seed") > 0) {
            return std::mt19937(vm["seed"].as<unsigned int>());
        } else {
            std::random_device rd;
            return std::mt19937(rd());
        }
    }();
    auto particles = initialParticles<ParticleT>(particlesNum, allSpecies, speciesInd, areaL, gen, particlesInit, particlesType);
    particlesNum = particles.size();

    auto force = Potential::force(vm);
    auto potential = Potential::potential(vm);

// ================================== output files ==================================

    std::string adiosOutput = vm["saveFile"].as<std::string>();

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("SimulationOutput");

    constexpr int D = degreesOfFreedom<ParticleT>();

    adios2::Variable<unsigned int> varStep = io.DefineVariable<unsigned int>("step");
    adios2::Variable<unsigned int> varVisitCounts = io.DefineVariable<unsigned int>("vistCount");
    adios2::Variable<unsigned int> varEndpointCounts = io.DefineVariable<unsigned int>("endpointCount");

    io.DefineAttribute<std::string>("particlesType", particlesType);
    io.DefineAttribute<Real>("particlesNum", particlesNum);
    io.DefineAttribute<Real>("areaL", areaL);
    
    adios2::Engine engine = io.Open(adiosOutput, adios2::Mode::Write);

    
    // ================================== Main Loop ==================================
    
    clock_t startTime = clock();

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
            handednessVec.clear();
            orientationVec.clear();
            positionsVec.clear();
            velocitiesVec.clear();
            for (auto &p : particles) {
                auto handedness = p.h;
                auto orientation = p.d;
                auto pos = getGeneralizedPositions(p);
                auto vel = getGeneralizedVelocities(p);

                handednessVec.push_back(handedness);
                orientationVec.push_back(orientation);
                
                for (int i = 0; i < D; ++i) {
                    positionsVec.push_back(pos[i]);
                    velocitiesVec.push_back(vel[i]);
                }
            }

            engine.Put(varHandedness, handednessVec.data());
            engine.Put(varOrientation, orientationVec.data());
            engine.Put(varPositions, positionsVec.data());
            engine.Put(varVelocities, velocitiesVec.data());
            writeStateStep = step + writeStateIntervalSteps;
        }
        if (step == writeEnergyStep) {
            auto kineticE = calKineticEnergy(particles, allSpecies);
            auto potentialE = calPotentialEnergy(particles, allSpecies, boxPBC, potential);

            auto neighborCount = calAveNeighborList(particles, neighborDistances, areaL);
            auto rotationalOrder = calRotationalOrder(particles);

            auto COMvelocity = calCOMvelocity(particles, allSpecies);
            auto COMangularVelocity = calAngularMomentum2D(particles, allSpecies, boxPBC);
            
            auto realTemperature = calInternalTemperature(particles, allSpecies);

            engine.Put(varKineticEnergy, kineticE.data());
            engine.Put(varPotentialEnergy, potentialE);

            engine.Put(varNeighborCount, neighborCount.data());
            engine.Put(varRotationalOrder, rotationalOrder);

            engine.Put(varComVelocity, COMvelocity.data());
            engine.Put(varComAngVelocity, COMangularVelocity);

            engine.Put(varRealTemperature, realTemperature.data());

            writeEnergyStep = step + writeEnergyIntervalSteps;
    
        }
        if constexpr (LETTUCE_THERMOSTAT != ThermostatID::None) {
            if (step == thermoStep) {
                thermostat(particles, allSpecies);


                thermoStep = step + thermoIntervalSteps;
                if constexpr (LETTUCE_THERMOSTAT == ThermostatID::VelocityScaling) {
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