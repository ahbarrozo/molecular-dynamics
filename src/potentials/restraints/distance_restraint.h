#pragma once

#include "potentials/potential.h"

class DistanceRestraintPotential : public Potential {
    float r_low_, r_up_;  // Flat-bottom potential
    float k_;
};