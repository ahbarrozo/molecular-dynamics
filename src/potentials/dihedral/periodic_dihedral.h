#pragma once

#include "potentials/potential.h"

class PeriodicDihedralPotential : public Potential {
    float k_;
    int n_;        // Periodicity
    float phi0_;   // Phase shift
};