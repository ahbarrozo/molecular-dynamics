#pragma once
#include <memory>
#include <vector>

#include "config.h"
#include "potentials/potential.h"
#include "vector3.h"
#include "helpers/neighbors_list.h"

class MDSystem {
    private:
        // Particle data (SoA layout)
        std::vector<Vector3> positions_;
        std::vector<Vector3> velocities_;
        std::vector<Vector3> forces_;
        std::map<std::pair<int, int>, std::unique_ptr<Potential>> potentials_;
        std::vector<real> masses_;
        std::vector<int> types_; // Atom types (for different elements)
        
        std::vector<Vector3> positions_at_last_update_;

        // System properties
        Vector3 box_size_;
        real neighbor_list_cutoff_;
        real skin_radius_ = static_cast<real>(0.3);
        real time_step_;
        size_t n_particles_;
        std::unique_ptr<NeightborsList> neighbors_list_;

        void build_neighbor_list();

        // Integration algorithms access only
        friend class VelocityVerlet;
        friend class Leapfrog;
        friend class Berendsen;
        friend class NoseHoover;
    
    public:
        // Constructor
        MDSystem();
        MDSystem(const Vector3& box_size, real time_step, real neighbor_list_cutoff_);

        void add_particle(Vector3 pos, Vector3 vel, real mass, int type);
        void apply_PBC();
        void compute_forces();
        void integrate();  // Update positions/velocities
        real kinetic_energy() const;
        bool needs_neighbors_list_update();
        real potential_energy() const;  
        real temperature() const;
        real total_energy() const;

        // Read-only access to system state
        const std::vector<Vector3>& positions() const;
        const std::vector<Vector3>& velocities() const;
        const std::vector<Vector3>& forces() const;
        const std::vector<real>& masses() const;
        size_t n_particles() const;
        const Vector3& box_size() const;
        real time_step() const;
};