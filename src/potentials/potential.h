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
        virtual void calculate(const Vector3& r_ij,
                            float r_mag,           // Distance magnitude
                            Vector3& force_i,
                            Vector3& force_j,
                            float& energy) const = 0;
        
        virtual real cutoff() const = 0;
        virtual std::string name() const = 0;
        virtual PotentialType type() const = 0;
};

