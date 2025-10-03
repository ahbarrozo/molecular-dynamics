#pragma once

#include "config.h"
#include "vector3.h"

enum class PotentialType {
    NON_BONDED,
    BOND,
    ANGLE,
    DIHEDRAL,
    IMPROPER,
    EXTERNAL
};

class Potential {
    public:
        virtual ~Potential() = default;     
        virtual void calculate(
                            real r2,           // Distance magnitude
                            real& force_mag,
                            real& energy) const = 0;
        
        virtual real cutoff() const = 0;
        virtual std::string name() const = 0;
        virtual PotentialType type() const = 0;
};

