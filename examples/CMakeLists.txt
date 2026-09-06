# =========================================
# Add executables
# =========================================
add_executable(chisurfmd_md mainMD.cpp)
add_executable(chisurfmd_aggregation mainMDaggregation.cpp)

# =========================================
# Link libraries
# =========================================
target_link_libraries(chisurfmd_md Boost::program_options adios2::cxx11)
target_link_libraries(chisurfmd_aggregation Boost::program_options adios2::cxx11)
