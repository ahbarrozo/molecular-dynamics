#include "leapfrog.h"

void Leapfrog::step(MDSystem& system) {
    const float dt = system.time_step();
    const float half_dt = 0.5f * dt;
    
    auto& positions = system.positions_;
    auto& velocities = system.velocities_;
    auto& forces = system.forces_;
    const auto& masses = system.masses_;
    
    if (first_step_) {
        // First step: need to initialize velocities at t - dt/2
        // v(-dt/2) = v(0) - a(0)*dt/2
        
        // First compute initial forces if not already computed
        system.compute_forces();
        
        // Shift velocities back by half timestep
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (size_t i = 0; i < system.n_particles(); ++i) {
            Vector3 accel = forces[i] / masses[i];
            velocities[i] -= accel * half_dt;
        }
        
        first_step_ = false;
    }
        
    // 1. Update velocities by full timestep: v(t+dt/2) = v(t-dt/2) + a(t)*dt
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < system.n_particles(); ++i) {
        Vector3 accel = forces[i] / masses[i];
        velocities[i] += accel * dt;  // Now at t + dt/2
    }
    
    // 2. Update positions: x(t+dt) = x(t) + v(t+dt/2)*dt
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < system.n_particles(); ++i) {
        positions[i] += velocities[i] * dt;
    }
    
    // 3. Apply periodic boundary conditions
    system.apply_PBC();
    
    // 4. Compute new forces at t+dt for next step
    system.compute_forces();
}

std::vector<Vector3> Leapfrog::get_synchronized_velocities(const MDSystem& system) const {
    const real half_dt = static_cast<real>(0.5) * system.time_step();
    const auto& velocities = system.velocities();
    const auto& forces = system.forces();
    const auto& masses = system.masses();
    
    std::vector<Vector3> v_sync(system.n_particles());
    
    // Interpolate: v(t) = v(t-dt/2) + a(t)*dt/2
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < system.n_particles(); ++i) {
        Vector3 accel = forces[i] / masses[i];
        v_sync[i] = velocities[i] + accel * half_dt;
    }
    
    return v_sync;
}