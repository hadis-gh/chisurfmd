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
#include "lettuce/LennardJones.h"
#include "lettuce/CircleDistribution.h"
#include "lettuce/CirclesIntersectionFuncs.h"

namespace po = boost::program_options;


int main(int argc, char* argv[]){

    po::variables_map vm;
    {
        po::options_description desc("Allowed Options");
        desc.add_options()
            ("help,h", "print help")
            ("particlesNum,n", po::value<int>()->default_value(10), "number of initial particles")
            ("time", po::value<double>()->default_value(10.), "max simulation")
            ("dt", po::value<double>()->default_value(.001), "integration step")
            ("stepT", po::value<double>()->default_value(1), "measurement interval")
            ("exclusionRadius", po::value<double>()->default_value(.8), "exclusion radius")
            ("areaL", po::value<double>()->default_value(20.), "simulation size")
            ("seed", po::value<unsigned int>(), "random seed")
            ("particleInit", po::value<std::string>()->default_value("DLA"), "particle intitialisation")

        ;

        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

        if(vm.count("help") > 0)
        {
            std::cout << desc << std::endl;
            return 0;
        }
    }

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
    const double dt = 0.001;
    const double cutoff = 20;
    const double Time = vm["time"].as<double>();
    const double areaL = vm["areaL"].as<double>();
    const double radius = vm["exclusionRadius"].as<double>();

    InitialParticlesConfiguration config = InitialParticlesConfiguration::DLA;
    if (vm["particleInit"].as<std::string>() == "RANDOM")
        config = InitialParticlesConfiguration::RANDOM_CIRCLES;
    IntegrationMethod integrationMethod = VELOCITY_VERLET; 

    LennardJonesForce<double> LJForce(epsilon, sigma, cutoff);


    std::vector<Particle<double>> particles = initialParticles (particlesNum, areaL, radius, gen, config);

    for (auto &p :particles){
        p.mass = mass; 
    }

    clock_t startTime = clock();

    const auto stepT = vm["stepT"].as<double>();
    int numSteps = static_cast<int>(Time / stepT);
    for (int a = 0; a < numSteps; ++a)
    {
        integrate(particles, dt, stepT, LJForce, integrationMethod);


    }

    clock_t endTime = clock();

    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    return 0;
}