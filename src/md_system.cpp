#include "md_system.h"

// Helper function to order pairs, used for storing potential functions
std::pair<std::string, std::string> ordered_pair(std::string a, std::string b) {
    return (a < b) ? 
           std::pair<std::string, std::string>{a, b} :
           std::pair<std::string, std::string>{b, a};
}

MDSystem::MDSystem(const Vector3& box_size, real time_step, real neighbor_list_cutoff)
    : box_size_(box_size),
      neighbor_list_cutoff_(neighbor_list_cutoff),
      time_step_(time_step) {

    n_particles_ = 0;
};

const std::vector<Vector3>& MDSystem::positions() const { return positions_; }
const std::vector<Vector3>& MDSystem::velocities() const { return velocities_; }
const std::vector<Vector3>& MDSystem::forces() const { return forces_; }
const real MDSystem::potential_energy() const { return potential_energy_; }
const std::vector<real>& MDSystem::masses() const { return masses_; }
const std::vector<std::string>& MDSystem::types() const { return types_; }
  
size_t MDSystem::n_particles() const { return n_particles_; }
const Vector3& MDSystem::box_size() const { return box_size_; }
real MDSystem::time_step() const { return time_step_; }

void MDSystem::add_particle(Vector3 pos, Vector3 vel, 
                            std::string type, real mass, 
                            real epsilon, real sigma) {
    positions_.push_back(pos);
    positions_at_last_update_.push_back(pos);
    velocities_.push_back(vel);
    masses_.push_back(mass);
    types_.push_back(type);
    forces_.push_back(Vector3(0.0f, 0.0f, 0.0f));
    
    for (auto& t: types_) {
        std::pair atom_pair = ordered_pair(type, t);
        if (potentials_.find(atom_pair) != potentials_.end()) continue;

        atom_types_parameters_[type] = {epsilon, sigma, mass};
        real epsilon_mixed = std::sqrt(atom_types_parameters_[type].epsilon * 
                                        atom_types_parameters_[t].epsilon);
        real sigma_mixed = static_cast<real>(0.5f) * (atom_types_parameters_[type].sigma + 
                                    atom_types_parameters_[t].sigma);
        
        potentials_[atom_pair] = std::make_unique<LennardJonesPotential>(
            epsilon_mixed, sigma_mixed, neighbor_list_cutoff_
        );
    }
    n_particles_++;
    atom_types_parameters_[type] = {epsilon, sigma, mass};
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
    real total_potential_energy = 0.0;

    if (!neighbors_list_) {
        neighbors_list_ = std::make_unique<NeightborsList>(box_size_, neighbor_list_cutoff_, n_particles_);
        neighbors_list_->update(positions_);
    }
    
    if (needs_neighbors_list_update()) {
        neighbors_list_->update(positions_);
    }
    
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < n_particles_; ++i) {
        forces_[i] = Vector3(0.0f, 0.0f, 0.0f);
    }

    // Compute forces directly using cell list
    #ifdef _OPENMP
        #pragma omp parallel reduction(+:total_potential_energy)
    #endif
    {
        std::vector<Vector3> thread_forces(n_particles_, Vector3(0, 0, 0));

        #ifdef _OPENMP
            #pragma omp for
        #endif
        for (size_t i = 0; i < n_particles_; ++i) {
    
            // Get neighbors on-the-fly
            auto neighbors = neighbors_list_->particles_in_cutoff(i, positions_, neighbor_list_cutoff_);
        
            for (size_t j : neighbors) {
                if (j > i) {
                    Vector3 dr = positions_[j] - positions_[i];
                    dr.wrap_PBC(box_size_);
                    real r2 = dr.mod2();
                    real force_mag, energy;

                    potentials_[ordered_pair(types_[i], types_[j])]->calculate(r2, force_mag, energy);
 
                    //compute_bonded_forces(); 
                    Vector3 force = dr * force_mag;
                
                    thread_forces[i] -= force;
                    thread_forces[j] += force;
                    total_potential_energy += energy;
                }
            }
        }
        
        // Combine thread-local forces into global forces
        #ifdef _OPENMP
        
        #pragma omp critical
        #endif
        {
            for (size_t i = 0; i < n_particles_; ++i) {
                forces_[i] += thread_forces[i];
            }
        }
    }
    
    potential_energy_ = total_potential_energy;
}

// Function to initialize velocities with Maxwell-Boltzmann distribution
void MDSystem::initialize_velocities(real initial_temperature) {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    for (size_t i = 0; i < n_particles_; ++i) {
        real sigma_v = std::sqrt(initial_temperature / masses_[i]);
        std::normal_distribution<float> dist(0.0f, sigma_v);
        
        velocities_[i] = Vector3(dist(gen), dist(gen), dist(gen));
    }
    
    // Remove center of mass motion
    Vector3 vcm(0.0f, 0.0f, 0.0f);
    float total_mass = 0.0f;
    
    for (size_t i = 0; i < n_particles_; ++i) {
        vcm += velocities_[i] * masses_[i];
        total_mass += masses_[i];
    }
    vcm = vcm / total_mass;
    
    for (size_t i = 0; i < n_particles_; ++i) {
        velocities_[i] -= vcm;
    }
}

real MDSystem::temperature() const {
    // Temperature from kinetic energy: T = 2*KE / (N_dof * kB)
    // In reduced units, kB = 1
    
    // Degrees of freedom = 3N - 3 (subtract COM motion)
    // For periodic boundaries, sometimes use 3N
    size_t n_dof = 3 * n_particles_ - 3;
    
    // If system is too small, avoid division issues
    if (n_particles_ < 2) {
        return 0.0f;
    }
    
    float ke = kinetic_energy();
    return 2.0f * ke / static_cast<float>(n_dof);
}

real MDSystem::kinetic_energy() const {
    float total_ke = 0.0f;
    
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:total_ke)
    #endif
    for (size_t i = 0; i < n_particles_; ++i) {
        real v2 = velocities_[i].mod2();
        total_ke += static_cast<real>(0.5f) * masses_[i] * v2;
    }
    
    return total_ke;
}

Vector3 MDSystem::center_of_mass_velocity() const {
    Vector3 p_total(0.0f, 0.0f, 0.0f);
    real m_total = 0.0f; 
    
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:m_total)
    #endif
    for (size_t i = 0; i < n_particles_; ++i) {
        m_total += masses_[i];
    }
    
    // Can't easily parallelize Vector3 reduction with OpenMP
    for (size_t i = 0; i < n_particles_; ++i) {
        p_total += velocities_[i] * masses_[i];
    }
    
    return p_total / m_total;
}

real MDSystem::kinetic_energy_center_of_mass() const {
    // Kinetic energy of the center of mass
    Vector3 v_com = center_of_mass_velocity();
    
    float m_total = 0.0f;
    for (size_t i = 0; i < n_particles_; ++i) {
        m_total += masses_[i];
    }
    
    return static_cast<real>(0.5f) * m_total * v_com.mod2();
}

real MDSystem::kinetic_energy_thermal() const {
    // Kinetic energy relative to center of mass (thermal motion only)
    Vector3 v_com = center_of_mass_velocity();
    real ke_thermal = 0.0f;
    
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:ke_thermal)
    #endif
    for (size_t i = 0; i < n_particles_; ++i) {
        Vector3 v_rel = velocities_[i] - v_com;
        ke_thermal += static_cast<real>(0.5f) * masses_[i] * v_rel.mod2();
    }
    
    return ke_thermal;
}

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
