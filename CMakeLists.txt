# ==============================================================================
# 1: PROJECT SETUP
# ==============================================================================
cmake_minimum_required(VERSION 3.20)
project(chisurfmd LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ==============================================================================
# 2: EXTERNAL DEPENDENCIES
# ==============================================================================
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Adios2 REQUIRED)

# ==============================================================================
# 3: USER CONFIGURATION (ccmake GUI Options)
# ==============================================================================
# 1. Particle Type Configuration
set(chisurfmd_PARTICLE_IDs "ParticleDot" "ParticleOriented")
set(chisurfmd_PARTICLE "ParticleOriented" CACHE STRING "Particle type")
set_property(CACHE chisurfmd_PARTICLE PROPERTY STRINGS ${chisurfmd_PARTICLE_IDs})

# 2. Thermostat Configuration
set(chisurfmd_THERMOSTAT_IDs "ThermostatID::VelocityScaling" "ThermostatID::Berendsen" "ThermostatID::Andersen" "ThermostatID::None")
set(chisurfmd_THERMOSTAT "ThermostatID::None" CACHE STRING "Thermostat")
set_property(CACHE chisurfmd_THERMOSTAT PROPERTY STRINGS ${chisurfmd_THERMOSTAT_IDs})

# 3. Potential Configuration
set(chisurfmd_POTENTIAL_IDs "IsotropicLJ" "OrientedLJ" "GeometricLJ" "TabularDFT" "PatchyLJ" "ChiralPatchyLJ" "ChiMorse")
set(chisurfmd_POTENTIAL "OrientedLJ" CACHE STRING "Potential")
set_property(CACHE chisurfmd_POTENTIAL PROPERTY STRINGS ${chisurfmd_POTENTIAL_IDs})

# ==============================================================================
# 4: COMPILER SETTINGS & MACROS
# ==============================================================================
add_definitions("-DCHISURFMD_PARTICLE=${chisurfmd_PARTICLE}")
add_definitions("-DCHISURFMD_THERMOSTAT=${chisurfmd_THERMOSTAT}")
add_definitions("-DCHISURFMD_POTENTIAL=${chisurfmd_POTENTIAL}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# ==============================================================================
# 5: DIRECTORIES & TARGETS
# ==============================================================================
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
)

add_subdirectory(examples)