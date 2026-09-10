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