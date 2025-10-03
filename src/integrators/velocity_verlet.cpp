#include "integrators/velocity_verlet.h"

void VelocityVerlet::step(MDSystem& system) {
    const float dt = system.time_step();
    const float dt2 = dt * dt;
    
    auto& positions = system.positions_;
    auto& velocities = system.velocities_;
    auto& forces = system.forces_;
    const auto& masses = system.masses();
    
    // Store old forces (needed for velocity update)
    std::vector<Vector3> old_forces = system.forces_;
    
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < system.n_particles(); ++i) {
        // Update positions: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
        Vector3 accel = forces[i] / masses[i];
        positions[i] += velocities[i] * dt + accel * (static_cast<real>(0.5) * dt2);
    }
    
    system.apply_PBC();
    system.compute_forces();
    
    // Update velocities: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t i = 0; i < system.n_particles(); ++i) {
        Vector3 accel_old = old_forces[i] / masses[i];
        Vector3 accel_new = forces[i] / masses[i];
        velocities[i] += (accel_old + accel_new) * (static_cast<real>(0.5) * dt);
    }
}