#pragma once

// Precision type for the entire MD simulation
#ifdef USE_DOUBLE_PRECISION
    using real = double;
    #define REAL_LITERAL(x) x
#else
    using real = float;
    #define REAL_LITERAL(x) x##f
#endif

// Optional: Other common type aliases
// using size_type = std::size_t;
// using index_type = std::size_t;
// using particle_id = std::size_t;

// Constants that depend on precision
namespace constants {
    constexpr real PI = static_cast<real>(3.14159265358979323846);
    constexpr real KB = static_cast<real>(1.380649e-23);  // J/K
}