#include <string>

#include "thermostats/berendsen.h"

Berendsen::Berendsen(real temp, real tau)
    : Thermostat(temp),
      tau_(tau) {
        if (tau <= 0) {
            throw std::invalid_argument("Tau must be positive");
        }
}

void Berendsen::apply(MDSystem& system) {
    real current_temp = system.temperature();
    
    // Avoid division by zero
    if (current_temp < 1e-10) {
        return;
    }
    
    // Berendsen scaling factor
    // λ = sqrt(1 + dt/tau * (T_target/T_current - 1))
    real dt = system.time_step();
    real scaling_factor = std::sqrt(1.0 + (dt / tau_) * 
                                      (target_temperature_ / current_temp - 1.0));
    
    // Apply scaling to all velocities
    auto& velocities = system.velocities_;
    
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < system.n_particles(); ++i) {
        velocities[i] *= scaling_factor;
    }
}