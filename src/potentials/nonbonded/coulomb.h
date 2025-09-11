#pragma once

#include "potentials/potential.h"

class CoulombPotential : public Potential {
    float dielectric_constant_;
    float cutoff_;
    // TODO: PME, Ewald, or reaction field
};