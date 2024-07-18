#pragma once

#include <boost/program_options.hpp>
#include <iostream>

//------------------------------------------------------------------
// Handle command line options:
//
// Options:
//  [-h|--help]                         Show help menue
//  [-n|--particlesNum]                 Number of initial particles
//  [-t|--time]                         Max simulation
//  [-dt]                               Integration step
//  [-writeInterval]                    Measurement interval
//  [-thermoInterval]                   Interval after which to apply thermostat
//  [-T|--temperature]                  Temperature (for the thermostat)
//  [-exclusionRadius]                  Exclusion radius
//  [-areaL]                            Simulation area size
//  [-seed]                             Random seed
//  [-particleInit]                     Particle initialization
//  [-]
//  [-]
//------------------------------------------------------------------

namespace po = boost::program_options;

po::variables_map parseCommandLineOptions(int argc, char* argv[]) {
    po::variables_map vm;
    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help,h", "print help")
        ("particlesNum,n", po::value<int>()->default_value(10), "number of initial particles")
        ("time,t", po::value<double>()->default_value(100.), "max simulation")
        ("dt", po::value<double>()->default_value(.001), "integration step")
        ("writeInterval", po::value<double>()->default_value(1.), "measurement interval")
        ("thermoInterval", po::value<double>()->default_value(1.), "interval after which to apply thermostat")
        ("temperature,T", po::value<double>()->default_value(1.), "temperature")
        ("exclusionRadius", po::value<double>()->default_value(.8), "exclusion radius")
        ("areaL", po::value<double>()->default_value(100.), "simulation size")
        ("seed", po::value<unsigned int>(), "random seed")
        ("particleInit", po::value<std::string>()->default_value("DLA"), "particle initialization");

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    if (vm.count("help") > 0) {
        std::cout << desc << std::endl;
        return 0;       //exit(0); ?
    }

    return vm;
}