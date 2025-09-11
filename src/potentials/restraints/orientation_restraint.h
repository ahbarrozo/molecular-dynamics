#pragma once

#include "potentials/potential.h"

class OrientationRestraintPotential : public Potential {
    Vector3 reference_vector_;
    float k_;
};