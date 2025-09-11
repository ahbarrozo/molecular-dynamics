#include <algorithm>

#include "helpers/neighbors_list.h"

NeightborsList::NeightborsList(const Vector3& box_size, real cutoff, size_t n_particles)
    : box_size_(box_size),
      cell_size_(cutoff),
      n_particles_(n_particles) {

    nx_ = static_cast<int>(box_size.x() / cell_size_);
    ny_ = static_cast<int>(box_size.y() / cell_size_);
    nz_ = static_cast<int>(box_size.z() / cell_size_);

    cell_size_ = std::min({
        box_size.x() / nx_, 
        box_size.y() / ny_, 
        box_size.z() / nz_
    });
    
    total_cells_ = nx_ * ny_ * nz_;
    cells_.resize(total_cells_);
}

int NeightborsList::get_cell_index(int ix, int iy, int iz) const {
    // Apply periodic boundary conditions to cell indices
    ix = (ix + nx_) % nx_;
    iy = (iy + ny_) % ny_;
    iz = (iz + nz_) % nz_;
    
    return ix + nx_ * (iy + ny_ * iz);
}

int NeightborsList::get_cell_index(const Vector3& pos) const {
    // Map position to cell indices
    // Handle positions outside box due to PBC
    real x = std::fmod(pos.x() + box_size_.x(), box_size_.x());
    real y = std::fmod(pos.y() + box_size_.y(), box_size_.y());
    real z = std::fmod(pos.z() + box_size_.z(), box_size_.z());
    
    int ix = static_cast<int>(x * nx_ / box_size_.x());
    int iy = static_cast<int>(y * ny_ / box_size_.y());
    int iz = static_cast<int>(z * nz_ / box_size_.z());
    
    // Clamp to valid range
    ix = std::min(nx_ - 1, std::max(0, ix));
    iy = std::min(ny_ - 1, std::max(0, iy));
    iz = std::min(nz_ - 1, std::max(0, iz));
    
    return get_cell_index(ix, iy, iz);
}

void NeightborsList::get_cell_coords(int cell_id, int& ix, int& iy, int& iz) const {
    iz = cell_id / (nx_ * ny_);
    iy = (cell_id - iz * nx_ * ny_) / nx_;
    ix = cell_id - iz * nx_ * ny_ - iy * nx_;
}

void NeightborsList::clear() {
    for (auto& cell : cells_) {
        cell.clear();
    }
}

void NeightborsList::update(const std::vector<Vector3>& positions) {
    
    // Estimate particles per cell
    size_t estimated_particles_per_cell = n_particles_ / total_cells_;

    // Clear all cells
    for (auto& cell : cells_) {
        cell.clear();
        cell.reserve(estimated_particles_per_cell);  // Avoid reallocations upon push_back
    }
    
    if (particle_cell_.size() != n_particles_) {
        particle_cell_.resize(n_particles_);
    }
    
    // Assign particles to cells
    for (size_t i = 0; i < n_particles_; ++i) {
        int cell_idx = get_cell_index(positions[i]);
        cells_[cell_idx].push_back(i);
        particle_cell_[i] = cell_idx;
    }
}

std::vector<size_t> NeightborsList::particles_in_neighbor_cells(size_t particle_idx) const {
    std::vector<size_t> neighbors;

    // Estimating the number of neighbours in average through all the eight 
    // adjacent cells, taking into account the cutoff stopping at half the cell
    size_t estimated_neighbors = (n_particles_ / total_cells_) * 13;

    neighbors.reserve(estimated_neighbors);
    
    // Get cell coordinates of particle
    int cell_idx = particle_cell_[particle_idx];
    int ix, iy, iz;
    get_cell_coords(cell_idx, ix, iy, iz);
    
    // Search 27 neighboring cells (including own cell)
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                int neighbor_cell = get_cell_index(ix + dx, iy + dy, iz + dz);
                
                // Add all particles in this cell
                for (size_t j : cells_[neighbor_cell]) {
                    if (j != particle_idx) {  // Don't include self
                        neighbors.push_back(j);
                    }
                }
            }
        }
    }
    
    return neighbors;
}

std::vector<size_t> NeightborsList::particles_in_cutoff(size_t particle_idx,
                                                  const std::vector<Vector3>& positions,
                                                  real cutoff) const {
    std::vector<size_t> neighbors;

    // This is a rough estimate
    size_t estimated_neighbors = (n_particles_ / total_cells_) * 2;

    neighbors.reserve(estimated_neighbors);
    
    const real cutoff2 = cutoff * cutoff;
    const Vector3& pos_i = positions[particle_idx];
    
    // Get cell coordinates
    int cell_idx = particle_cell_[particle_idx];
    int ix, iy, iz;
    get_cell_coords(cell_idx, ix, iy, iz);
    
    // Search neighboring cells
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                int neighbor_cell = get_cell_index(ix + dx, iy + dy, iz + dz);
                
                // Check all particles in this cell
                for (size_t j : cells_[neighbor_cell]) {
                    if (j != particle_idx) {
                        Vector3 dr = positions[j] - pos_i;
                        dr.wrap_PBC(box_size_);
                        
                        if (dr.mod2() < cutoff2) {
                            neighbors.push_back(j);
                        }
                    }
                }
            }
        }
    }
    
    return neighbors;
}