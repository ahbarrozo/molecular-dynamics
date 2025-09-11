#pragma once

#include "potentials/potential.h"

class ImproperDihedralPotential : public Potential {
    float k_;
    float xi0_;  // Equilibrium improper angle
};