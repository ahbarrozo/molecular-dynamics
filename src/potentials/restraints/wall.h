#pragma once

#include "potentials/potential.h"

class WallPotential : public Potential {
    enum WallType { HARMONIC, LJ_9_3, EXPONENTIAL };
    Vector3 normal_;
    float position_;
    float k_;
};