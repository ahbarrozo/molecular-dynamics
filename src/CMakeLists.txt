add_library(md_lib
    vector3.cpp
    vector3.h
    
    helpers/neighbors_list.cpp
    helpers/neighbors_list.h
    integrators/integrator.h
    integrators/leapfrog.cpp
    integrators/leapfrog.h
    integrators/velocity_verlet.cpp
    integrators/velocity_verlet.h
    potentials/potential.h
    potentials/angle/harmonic_angle.h
    potentials/angle/harmonic_cosine_angle.h
    potentials/angle/urey_bradley.h
    potentials/bond/harmonic_bond.h
    potentials/bond/morse.h
    potentials/bond/quartic_bond.h
    potentials/dihedral/improper_dihedral.h
    potentials/dihedral/opls_dihedral.h
    potentials/dihedral/periodic_dihedral.h
    potentials/dihedral/ryckaert_bellemans.h
    potentials/nonbonded/buckingham.h
    potentials/nonbonded/coulomb.h
    potentials/nonbonded/lennard_jones.h
    potentials/nonbonded/soft_core.h
    potentials/restraints/distance_restraint.h
    potentials/restraints/orientation_restraint.h
    potentials/restraints/position_restraint.h
    potentials/restraints/wall.h
    thermostats/thermostat.h
    thermostats/berendsen.cpp
    thermostats/berendsen.h
    thermostats/nose_hoover.cpp
    thermostats/nose_hoover.h
    md_system.cpp
    md_system.h
)

# Link OpenBLAS to the library
target_link_libraries(md_lib PUBLIC ${OPENBLAS_LIBRARIES})
target_include_directories(md_lib PUBLIC ${OPENBLAS_INCLUDE_DIRS})
target_compile_options(md_lib PUBLIC ${OPENBLAS_CFLAGS_OTHER})

# Make headers available to other parts of project and external users
target_include_directories(md_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Set C++ standard for the library
target_compile_features(md_lib PUBLIC cxx_std_17)

# Add compiler warnings for better code quality
target_compile_options(md_lib PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

option(USE_DOUBLE_PRECISION "Use double precision instead of float" OFF)

if(USE_DOUBLE_PRECISION)
    target_compile_definitions(md_lib PUBLIC USE_DOUBLE_PRECISION)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(md_lib PUBLIC OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found - parallel execution enabled")
else()
    message(STATUS "OpenMP not found - using sequential execution")
endif()

# Link the executable with your library (no need for separate tensor_lib)
# target_link_libraries(train_network PRIVATE md_lib)

# Optional: Add a test executable if you have tests
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_vector3.cpp)
    add_executable(test_vector3 test_vector3.cpp)
    target_link_libraries(test_vector3 PRIVATE md_lib)
endif()