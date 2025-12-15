cmake_minimum_required(VERSION 3.20)
# enable_testing()
project(lettuce LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

#find_package(Boost CONFIG)
find_package(Boost COMPONENTS program_options REQUIRED)
set(Adios2_DIR "/opt/pkg/adios2/lib/cmake/adios2")
find_package(Adios2 REQUIRED)

set(lettuce_md_PARTICLE_IDs
	"ParticleDot"
	"ParticleOriented"
	)
set(lettuce_md_PARTICLE "ParticleOriented" CACHE STRING "Particle type")
set_property(CACHE lettuce_md_PARTICLE PROPERTY STRINGS ${lettuce_md_PARTICLE_IDs})

set(lettuce_md_THERMOSTAT_IDs
	"ThermostatID::VelocityScaling"
	"ThermostatID::Berendsen"
	"ThermostatID::Andersen"
	"ThermostatID::None"
	)
set(lettuce_md_THERMOSTAT "ThermostatID::None" CACHE STRING "Thermostat")
set_property(CACHE lettuce_md_THERMOSTAT PROPERTY STRINGS ${lettuce_md_THERMOSTAT_IDs})

set(lettuce_POTENTIAL_IDs
	"IsotropicLJ"
	"OrientedLJ"
	"GeometricLJ"
	"TabularDFT"
	"PatchyLJ"
	"ChiralPatchyLJ"
	)
set(lettuce_POTENTIAL "IsotropicLJ" CACHE STRING "Potential")
set_property(CACHE lettuce_POTENTIAL PROPERTY STRINGS ${lettuce_POTENTIAL_IDs})

add_definitions("-DLETTUCE_PARTICLE=${lettuce_md_PARTICLE}")
add_definitions("-DLETTUCE_THERMOSTAT=${lettuce_md_THERMOSTAT}")
add_definitions("-DLETTUCE_POTENTIAL=${lettuce_POTENTIAL}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${Boost_INCLUDE_DIRS}
	#${Adios2_DIR}/include
)

# add_subdirectory(thirdParty/lettuce)

#add_subdirectory(include)
#add_subdirectory(src)

option(BUILD_TESTING "Build the testing tree." OFF)
enable_testing()
if(BUILD_TESTING)
  add_subdirectory(test/unit)	
else()
	add_subdirectory(test)
endif()
