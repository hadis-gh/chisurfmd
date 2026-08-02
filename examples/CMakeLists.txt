# =========================================
# Add executables
# =========================================
add_executable(testMD mainMD.cpp)
add_executable(testMDaggregation mainMDaggregation.cpp)

# =========================================
# Link libraries
# =========================================
target_link_libraries(testMD Boost::program_options adios2::cxx11)
target_link_libraries(testMDaggregation Boost::program_options adios2::cxx11)
