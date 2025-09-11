#pragma once

#include "potentials/potential.h"

class RyckaertBellemansPotential : public Potential {
    float C0_, C1_, C2_, C3_, C4_, C5_;  // Coefficients
};