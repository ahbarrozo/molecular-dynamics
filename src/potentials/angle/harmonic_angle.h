#pragma once

#include "potentials/potential.h"

class HarmonicAnglePotential : public Potential {
    float k_;
    float theta0_;  // Equilibrium angle
};