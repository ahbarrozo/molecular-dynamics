#include "nose_hoover.h"
#include <cmath>

// Simple Nosé-Hoover
NoseHooverThermostat::NoseHooverThermostat(real temp, real tau)
    : Thermostat(temp), 
      tau_(tau),
      xi_(0.0),
      xi_dot_(0.0),
      use_chains_(false),
      chain_length_(0) {
    
    if (tau <= 0) {
        throw std::invalid_argument("Tau must be positive");
    }
}

// Nosé-Hoover chains
NoseHooverThermostat::NoseHooverThermostat(real temp, real tau, int chain_length)
    : Thermostat(temp),
      tau_(tau),
      xi_(0.0),
      xi_dot_(0.0),
      use_chains_(true),
      chain_length_(chain_length) {
    
    if (tau <= 0) {
        throw std::invalid_argument("Tau must be positive");
    }
    if (chain_length < 2) {
        throw std::invalid_argument("Chain length must be at least 2");
    }
    
    xi_chain_.resize(chain_length, 0.0);
    xi_dot_chain_.resize(chain_length, 0.0);
}

void NoseHooverThermostat::apply(MDSystem& system) {
    const real dt = system.time_step();
    const real half_dt = static_cast<real>(0.5) * dt;
    const size_t n_particles = system.n_particles();
    const size_t n_dof = 3 * n_particles - 3;  // Degrees of freedom (subtract 3 for COM motion)
    
    auto& velocities = system.velocities_;
    const auto& masses = system.masses();
    
    // Calculate thermal mass Q = n_dof * kB * T * tau^2
    // In reduced units where kB = 1:
    Q_ = n_dof * target_temperature_ * tau_ * tau_;
    
    if (!use_chains_) {
        // Simple Nosé-Hoover thermostat
        
        // Calculate current kinetic energy
        real kinetic_energy = static_cast<real>(0.0);
        #ifdef _OPENMP
            #pragma omp parallel for reduction(+:kinetic_energy)
        #endif
        for (size_t i = 0; i < n_particles; ++i) {
            real v2 = velocities[i].mod2();
            kinetic_energy += static_cast<real>(0.5) * masses[i] * v2;
        }
        
        // Update xi_dot (acceleration of the thermostat variable)
        // xi_dot = (2*KE - n_dof*kB*T) / Q
        xi_dot_ = (static_cast<real>(2.0) * kinetic_energy - 
                   static_cast<real>(n_dof) * target_temperature_) / Q_;
        
        // Update xi (half step)
        xi_ += xi_dot_ * half_dt;
        
        // Scale velocities: v = v * exp(-xi * dt)
        real scaling_factor = std::exp(-xi_ * dt);
        
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (size_t i = 0; i < n_particles; ++i) {
            velocities[i] *= scaling_factor;
        }
        
        // Recalculate kinetic energy after scaling
        kinetic_energy *= scaling_factor * scaling_factor;
        
        // Update xi_dot with new kinetic energy
        xi_dot_ = (2.0f * kinetic_energy - 
                   n_dof * target_temperature_) / Q_;
        
        // Complete xi update (second half step)
        xi_ += xi_dot_ * half_dt;
        
    } else {
        // Nosé-Hoover chains (better temperature control)
        
        // Calculate current kinetic energy
        real kinetic_energy = static_cast<real>(0.0);
        #ifdef _OPENMP
            #pragma omp parallel for reduction(+:kinetic_energy)
        #endif
        for (size_t i = 0; i < n_particles; ++i) {
            real v2 = velocities[i].mod2();
            kinetic_energy += static_cast<real>(0.5) * masses[i] * v2;
        }
        
        // Update chain variables (from end to beginning)
        for (int j = chain_length_ - 1; j >= 0; --j) {
            real Q_j = (j == 0) ? Q_ : target_temperature_;
            
            if (j == 0) {
                // First chain element couples to the system
                xi_dot_chain_[0] = (static_cast<real>(2.0) * kinetic_energy - 
                                   static_cast<real>(n_dof) * target_temperature_) / Q_j;
            } else if (j == chain_length_ - 1) {
                // Last chain element
                xi_dot_chain_[j] = (xi_chain_[j-1] * xi_chain_[j-1] - target_temperature_) / Q_j;
            } else {
                // Middle chain elements
                xi_dot_chain_[j] = (xi_chain_[j-1] * xi_chain_[j-1] - target_temperature_) / Q_j;
                xi_dot_chain_[j] -= xi_chain_[j+1] * xi_dot_chain_[j];
            }
            
            // Update xi for this chain element
            xi_chain_[j] += xi_dot_chain_[j] * half_dt;
        }
        
        // Scale velocities using first chain element
        real scaling_factor = std::exp(-xi_chain_[0] * dt);
        
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (size_t i = 0; i < n_particles; ++i) {
            velocities[i] *= scaling_factor;
        }
        
        // Complete chain updates (second half step)
        kinetic_energy *= scaling_factor * scaling_factor;
        
        for (int j = 0; j < chain_length_; ++j) {
            real Q_j = (j == 0) ? Q_ : target_temperature_;
            
            if (j == 0) {
                xi_dot_chain_[0] = (static_cast<real>(2.0) * kinetic_energy - 
                                   static_cast<real>(n_dof) * target_temperature_) / Q_j;
            }
            
            xi_chain_[j] += xi_dot_chain_[j] * half_dt;
        }
        
        // Store first chain element in xi_ for monitoring
        xi_ = xi_chain_[0];
    }
}