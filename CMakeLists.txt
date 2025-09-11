cmake_minimum_required(VERSION 3.16)
project(MolecularDynamics VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Find OpenBLAS
find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENBLAS REQUIRED openblas)

# Find OpenMP
find_package(OpenMP)


# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)