add_executable(
    test_core
    test_core.cpp
)

target_include_directories(
    test_core
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

add_test(
    NAME CoreUtilities
    COMMAND test_core
)

#  ============================================
add_executable(
    test_oriented_lj
    test_oriented_lj.cpp
)

target_include_directories(
    test_oriented_lj
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

add_test(
    NAME OrientedLJ
    COMMAND test_oriented_lj
)

#  ============================================
add_executable(
    test_integration
    test_integration.cpp
)

target_include_directories(
    test_integration
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

add_test(
    NAME VelocityVerletNVE
    COMMAND test_integration
)