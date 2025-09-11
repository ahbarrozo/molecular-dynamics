#pragma once

#include "potentials/potential.h"

class HarmonicBondPotential : public Potential {
    float k_;   // Spring constant
    float r0_;  // Equilibrium distance
};