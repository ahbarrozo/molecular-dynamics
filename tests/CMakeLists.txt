# Create test executable
add_executable(test_vector3 test_vector3.cpp)

# Link with your tensor library (automatically gets OpenBLAS too)
target_link_libraries(test_vector3 PRIVATE md_lib)