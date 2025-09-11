#pragma once

#include "potentials/potential.h"

class QuarticBondPotential : public Potential {
    float k2_, k3_, k4_;  // Higher order terms
    float r0_;
};