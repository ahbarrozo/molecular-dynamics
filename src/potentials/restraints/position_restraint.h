#pragma once

#include "potentials/potential.h"

class PositionRestraintPotential : public Potential {
    Vector3 reference_position_;
    float k_x_, k_y_, k_z_;
};