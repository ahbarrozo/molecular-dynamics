#pragma once

#include "potentials/potential.h"
#include "config.h"

class LennardJonesPotential : public Potential {
    private:
        real epsilon_, sigma_, cutoff_, cutoff2_;
        real sigma6_;
        real energy_cutoff_;

    public: 
        LennardJonesPotential(real epsilon, real sigma, real cutoff);
        void calculate(real r2, real& force_mag, real& energy) const override;
        real cutoff() const override;
        std::string name() const override { return "Lennard-Jones"; }
        PotentialType type() const override { return PotentialType::NON_BONDED; }
};