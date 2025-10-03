#pragma once
#include <map>
#include <memory>
#include <random>
#include <vector>

#include "config.h"
#include "helpers/neighbors_list.h"
#include "potentials/potential.h"
#include "potentials/nonbonded/lennard_jones.h"
#include "vector3.h"

struct AtomTypeParameters {
    real epsilon;
    real sigma;
    real mass;
};

class MDSystem {
    private:
        // Particle data (SoA layout)
        std::vector<Vector3> positions_;
        std::vector<Vector3> velocities_;
        std::vector<Vector3> forces_;
        std::vector<real> masses_;
        std::vector<std::string> types_;
        std::map<std::pair<std::string, std::string>, std::unique_ptr<Potential>> potentials_;
        std::map<std::string, AtomTypeParameters> atom_types_parameters_;
        
        std::vector<Vector3> positions_at_last_update_;

        // System properties
        Vector3 box_size_;
        real neighbor_list_cutoff_;
        real skin_radius_ = static_cast<real>(0.3);
        real time_step_;
        size_t n_particles_;
        std::unique_ptr<NeightborsList> neighbors_list_;
        real potential_energy_;

        void build_neighbor_list();

        // std::unique_ptr<Integrator> integrator_;
        // std::unique_ptr<Thermostat> thermostat_;

        // Integration algorithms access only
        friend class VelocityVerlet;
        friend class Leapfrog;
        friend class BerendsenThermostat;
        friend class NoseHooverThermostat;
    
    public:
        // Constructor
        MDSystem();
        MDSystem(const Vector3& box_size, 
                 real time_step, 
                 real neighbor_list_cutoff_);

        void add_particle(Vector3 pos, Vector3 vel, 
                          std::string type, real mass, 
                          real epsilon, real sigma);
    
        void apply_PBC();
        Vector3 center_of_mass_velocity() const;
        void compute_forces();
        void initialize_velocities(real initial_temperature);    
        real kinetic_energy() const;
        real kinetic_energy_center_of_mass() const;
        real kinetic_energy_thermal() const;
        bool needs_neighbors_list_update();
        real temperature() const;
        real total_energy() const;
        
        // Read-only access to system state
        const std::vector<Vector3>& positions() const;
        const std::vector<Vector3>& velocities() const;
        const std::vector<Vector3>& forces() const;
        const std::vector<real>& masses() const;
        const std::vector<std::string>& types() const;
        const real potential_energy() const;  
        size_t n_particles() const;
        const Vector3& box_size() const;
        real time_step() const;
};