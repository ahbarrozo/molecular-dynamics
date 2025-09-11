#pragma once

#include "potentials/potential.h"

class SoftCorePotential : public Potential {
    float lambda_;  // Coupling parameter
    float alpha_;   // Softness parameter
};