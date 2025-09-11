#pragma once

#include "potentials/potential.h"

class OPLSDihedralPotential : public Potential {
    float k1_, k2_, k3_, k4_;  // Fourier coefficients
};