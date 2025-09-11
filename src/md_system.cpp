#include "md_system.h"

// Helper function to order pairs, used for storing potential functions
std::pair<int, int> ordered_pair(int a, int b) {
    return (a < b) ? std::pair<int, int>{a, b} : std::pair<int, int>{b, a};
}

MDSystem::MDSystem(const Vector3& box_size, real time_step, real neighbor_list_cutoff)
    : box_size_(box_size),
      neighbor_list_cutoff_(neighbor_list_cutoff),
      time_step_(time_step) {

    // TODO : Initialize system
};

const std::vector<Vector3>& MDSystem::positions() const { return positions_; }
const std::vector<Vector3>& MDSystem::velocities() const { return velocities_; }
const std::vector<Vector3>& MDSystem::forces() const { return forces_; }
const std::vector<float>& MDSystem::masses() const { return masses_; }
  
size_t MDSystem::n_particles() const { return n_particles_; }
const Vector3& MDSystem::box_size() const { return box_size_; }
float MDSystem::time_step() const { return time_step_; }

void MDSystem::add_particle(Vector3 pos, Vector3 vel, real mass, int type) {
    positions_.push_back(pos);
    positions_at_last_update_.push_back(pos);
    velocities_.push_back(vel);
    masses_.push_back(mass);
    forces_.push_back(Vector3(0.0f, 0.0f, 0.0f));
    
    for (auto& t: types_) {
        if (potentials_.find(ordered_pair(type, t)) {
            continue;
        }

        potentials_[ordered_pair(type, t)] = ...;
    }
    n_particles_++;
    types_.push_back(type);
}

void MDSystem::apply_PBC() {
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < n_particles_; ++i) {
        positions_[i].wrap_PBC(box_size_);
    }
}

void MDSystem::compute_forces() {
    if (!neighbors_list_) {
        neighbors_list_ = std::make_unique<NeightborsList>(box_size_, neighbor_list_cutoff_, n_particles_);
    }
    
    // Update if needed
    if (needs_neighbors_list_update()) {
        neighbors_list_->update(positions_);
    }
    
    // Compute forces directly using cell list
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < n_particles_; ++i) {
        forces_[i] = Vector3(0.0f, 0.0f, 0.0f);
        real total_potential_energy = 0.0;

    
        // Get neighbors on-the-fly
        auto neighbors = neighbors_list_->particles_in_cutoff(i, positions_, neighbor_list_cutoff_);
        
        for (size_t j : neighbors) {
            if (j > i) {
                Vector3 dr = positions_[j] - positions_[i];
                dr.wrap_PBC(box_size_);
                real r2 = dr.mod2();
            
                real force_mag, energy;
            
                // Calculate based on atom types
                if (mixed_potential_) {
                    mixed_potential_->calculate_for_types(types_[i], types_[j], dr, r2, force_mag, energy);
                } else if (pair_potential_) {
                    pair_potential_->calculate(dr, r2, force_mag, energy);
                }
            
                Vector3 force = dr * force_mag;
            
                // Ensure thread-safe updates in force increments
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                forces_[i].x() -= force.x();
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                forces_[i].y() -= force.y();
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                forces_[i].z() -= force.z();
            
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                forces_[j].x() += force.x();
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                forces_[j].y() += force.y();
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                forces_[j].z() += force.z();
            
            total_potential_energy += energy;
        }
    }
    
    potential_energy_ = total_potential_energy;
            }
        }
    }
}


    
    
        
        for (size_t j : neighbors) {
            if (j <= i) continue;  // Avoid double counting
            
 
    
    // Add bonded forces (if any)
    compute_bonded_forces();

// Checks whether one of the particles moved too much since the last update
bool MDSystem::needs_neighbors_list_update() {    
    const real threshold = skin_radius_ * static_cast<real>(0.5);
    const real threshold2 = threshold * threshold;
    
    // Check max displacement on demand
    for (size_t i = 0; i < n_particles_; ++i) {
        Vector3 disp = positions_[i] - positions_at_last_update_[i];
        if (disp.mod2() > threshold2) {
            positions_at_last_update_ = positions_;
            return true;
        }
    }
    return false;
}


