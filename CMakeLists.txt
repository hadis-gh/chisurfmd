# ==============================================================================
# 1: PROJECT SETUP
# ==============================================================================
cmake_minimum_required(VERSION 3.20)
project(lettuce LANGUAGES CXX)
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

set(Adios2_DIR "/opt/pkg/adios2/lib/cmake/adios2")
find_package(Adios2 REQUIRED)

# ==============================================================================
# 3: USER CONFIGURATION (ccmake GUI Options)
# ==============================================================================
# 1. Particle Type Configuration
set(lettuce_md_PARTICLE_IDs "ParticleDot" "ParticleOriented")
set(lettuce_md_PARTICLE "ParticleOriented" CACHE STRING "Particle type")
set_property(CACHE lettuce_md_PARTICLE PROPERTY STRINGS ${lettuce_md_PARTICLE_IDs})

# 2. Thermostat Configuration
set(lettuce_md_THERMOSTAT_IDs "ThermostatID::VelocityScaling" "ThermostatID::Berendsen" "ThermostatID::Andersen" "ThermostatID::None")
set(lettuce_md_THERMOSTAT "ThermostatID::None" CACHE STRING "Thermostat")
set_property(CACHE lettuce_md_THERMOSTAT PROPERTY STRINGS ${lettuce_md_THERMOSTAT_IDs})

# 3. Potential Configuration
set(lettuce_POTENTIAL_IDs "IsotropicLJ" "OrientedLJ" "GeometricLJ" "TabularDFT" "PatchyLJ" "ChiralPatchyLJ" "ChiMorse")
set(lettuce_POTENTIAL "OrientedLJ" CACHE STRING "Potential")
set_property(CACHE lettuce_POTENTIAL PROPERTY STRINGS ${lettuce_POTENTIAL_IDs})

# ==============================================================================
# 4: COMPILER SETTINGS & MACROS
# ==============================================================================
add_definitions("-DLETTUCE_PARTICLE=${lettuce_md_PARTICLE}")
add_definitions("-DLETTUCE_THERMOSTAT=${lettuce_md_THERMOSTAT}")
add_definitions("-DLETTUCE_POTENTIAL=${lettuce_POTENTIAL}")

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