#pragma once

#include "potentials/potential.h"

class CosineAnglePotential : public Potential {
    float k_;
    float cos_theta0_;
};