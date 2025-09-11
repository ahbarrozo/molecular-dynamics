#pragma once

#include "potentials/potential.h"

class BuckinghamPotential : public Potential {
    float A_, B_, C_;  // V = A*exp(-B*r) - C/r^6
};