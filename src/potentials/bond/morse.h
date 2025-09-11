#pragma once

#include "potentials/potential.h"

class MorseBondPotential : public Potential {
    float D_e_, a_, r_e_;
};
