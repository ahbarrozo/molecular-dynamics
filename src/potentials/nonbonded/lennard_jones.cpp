#include "potentials/nonbonded/lennard_jones.h"

LennardJonesPotential::LennardJonesPotential(real epsilon, real sigma, real cutoff)
    : epsilon_(epsilon), sigma_(sigma), cutoff_(cutoff) {
    
    cutoff2_ = cutoff * cutoff;
    sigma6_ = sigma * sigma * sigma * sigma * sigma * sigma;
    
    // Calculate energy shift for continuity at cutoff
    real rc2_inv = static_cast<real>(1.0) / cutoff2_;
    real rc6_inv = rc2_inv * rc2_inv * rc2_inv;
    real sigma6_rc6 = sigma6_ * rc6_inv;
    energy_cutoff_ = static_cast<real>(4.0) * epsilon_ * (sigma6_rc6 * sigma6_rc6 - sigma6_rc6);
}

void LennardJonesPotential::calculate(real r2, 
    real& force_mag, real& energy) const {

        if (r2 > cutoff2_) {
            force_mag = 0.0f;
            energy = 0.0f;
            return;
        }
        
        real r2_inv = static_cast<real>(1.0) / r2;
        real r6_inv = r2_inv * r2_inv * r2_inv;
        real sigma6_r6 = sigma6_ * r6_inv;
        real sigma12_r12 = sigma6_r6 * sigma6_r6;
        
        // Force magnitude divided by r
        force_mag = static_cast<real>(24.0) * epsilon_ * r2_inv * 
                   (static_cast<real>(2.0) * sigma12_r12 - sigma6_r6);
        
        // Shifted potential energy
        energy = static_cast<real>(4.0) * epsilon_ * (sigma12_r12 - sigma6_r6) - energy_cutoff_;
    }
    
real LennardJonesPotential::cutoff() const { return cutoff_; }