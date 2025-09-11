#pragma once

#include "potentials/potential.h"

class UreyBradleyPotential : public Potential {
    float k_angle_, theta0_;
    float k_ub_, r_ub_;  // 1-3 distance term
};