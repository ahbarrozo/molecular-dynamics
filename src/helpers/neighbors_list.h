#pragma once

#include <vector>

#include "config.h"
#include "vector3.h"

class NeightborsList {
    private:
        Vector3 box_size_;
        real cell_size_;
        size_t n_particles_;

        // Grid dimensions
        int nx_, ny_, nz_, total_cells_;

        std::vector<std::vector<size_t>> cells_; // Particle indexes in each cell
        std::vector<size_t> particle_cell_;  // Which cell each particle is in

        // Helper functions
        int get_cell_index(int ix, int iy, int iz) const;
        int get_cell_index(const Vector3& pos) const;
        void get_cell_coords(int cell_id, int& ix, int& iy, int& iz) const;
        
    public:
        NeightborsList(const Vector3& box_size, real cutoff, size_t n_particles);
        void update(const std::vector<Vector3>& positions);
        std::vector<size_t> particles_in_cutoff(size_t particle_id, 
                                                    const std::vector<Vector3>& positions,
                                                    real cutoff) const;
        std::vector<size_t> particles_in_neighbor_cells(size_t particle_id) const;

        void clear();
};