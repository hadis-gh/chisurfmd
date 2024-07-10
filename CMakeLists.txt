cmake_minimum_required(VERSION 3.20)
# enable_testing()
project(lettuce LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

#find_package(Boost CONFIG)
find_package(Boost COMPONENTS program_options REQUIRED)
#find_package(Adios2)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror")
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${Boost_INCLUDE_DIRS}
)
#add_subdirectory(include)
#add_subdirectory(src)
add_subdirectory(test)