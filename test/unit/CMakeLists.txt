
add_library(Catch2Main include/catch_amalgamated.cpp )
target_include_directories(Catch2Main PUBLIC include)


add_executable(test_initParticles testInitParticles.cpp)
add_executable(test_vec2D testVec2D.cpp)
add_executable(test_particleDot testParticleDot.cpp)

target_link_libraries(test_initParticles adios2::cxx11)
target_link_libraries(test_particleDot PRIVATE Catch2Main)


add_test(InitParticles test_initParticles)
add_test(VecClass2D test_vec2D)
add_test(ParticleDot test_particleDot)