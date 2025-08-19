
add_executable(test_initParticles testInitParticles.cpp)
add_executable(test_vec2D testVec2D.cpp)

target_link_libraries(test_initParticles adios2::cxx11)


add_test(InitParticles test_initParticles)
add_test(VecClass2D test_vec)