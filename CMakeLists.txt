cmake_minimum_required(VERSION 4.0.3)
project(VFHETest)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.11.0
)
FetchContent_MakeAvailable(googletest)

add_executable(vfhe_test vfhe_test.cpp)
target_link_libraries(vfhe_test gtest gtest_main)