// Test with 100 noble gas particles

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

#include "vector3.h"
#include "md_system.h"
#include "integrators/velocity_verlet.h"
#include "thermostats/berendsen.h"
#include "potentials/nonbonded/lennard_jones.h"

void compute_rdf(const MDSystem& system, const std::string& filename) {
    const int n_bins = 100;
    const real max_r = system.box_size().x() * 0.5f;  // Half box size
    const real dr = max_r / n_bins;
    
    std::vector<real> g_r(n_bins, 0.0f);
    const auto& positions = system.positions();
    const auto& types = system.types();
    const auto& box_size = system.box_size();
    
    for (size_t i = 0; i < system.n_particles(); ++i) {
        for (size_t j = i + 1; j < system.n_particles(); ++j) {
            Vector3 dr_vec = positions[j] - positions[i];
            dr_vec.wrap_PBC(box_size);
            real r = dr_vec.mod();
            
            if (r < max_r) {
                int bin = static_cast<int>(r / dr);
                g_r[bin] += 2.0f;  // Count for both i-j and j-i
            }
        }
    }
    
    // Normalize
    float n_density = system.n_particles() / 
                     (box_size.x() * box_size.y() * box_size.z());
    
    std::ofstream rdf_file(filename);
    rdf_file << "# r  g(r)\n";
    
    for (int i = 0; i < n_bins; ++i) {
        real r = static_cast<real>(i + 0.5f) * dr;
        real shell_volume = 4.0f * M_PI * r * r * dr;
        real expected_pairs = shell_volume * n_density * system.n_particles();
        g_r[i] /= expected_pairs;
        
        rdf_file << r << " " << g_r[i] << "\n";
    }
    rdf_file.close();
}

// Function to create initial positions on a lattice
void create_fcc_lattice(MDSystem& system, int n_cells, float lattice_const) {
    // FCC unit cell positions
    std::vector<Vector3> unit_cell = {
        Vector3(0.0f, 0.0f, 0.0f),
        Vector3(0.5f, 0.5f, 0.0f),
        Vector3(0.5f, 0.0f, 0.5f),
        Vector3(0.0f, 0.5f, 0.5f)
    };
    
    int particle_count = 0;
    for (int ix = 0; ix < n_cells; ++ix) {
        for (int iy = 0; iy < n_cells; ++iy) {
            for (int iz = 0; iz < n_cells; ++iz) {
                for (const auto& base_pos : unit_cell) {
                    if (particle_count >= 100) return;
                    
                    Vector3 pos = Vector3(
                        (ix + base_pos.x()) * lattice_const,
                        (iy + base_pos.y()) * lattice_const,
                        (iz + base_pos.z()) * lattice_const
                    );
                    
                    if (particle_count % 5 == 0)
                        system.add_particle(
                            pos, Vector3(0, 0, 0), 
                            "Ar", static_cast<real>(39.95f), 
                            1.0f, 1.0f);
                    else
                        system.add_particle(
                            pos, Vector3(0, 0, 0), 
                            "Ne", static_cast<real>(20.18f), 
                            static_cast<real>(0.5f), static_cast<real>(0.85f));

                    particle_count++;
                }
            }
        }
    }
}

int main() {
    // Simulation parameters
    const real dt = static_cast<real>(0.001f);           // Time step (reduced units)
    const real temperature = static_cast<real>(1.0f);    // Target temperature (reduced units)
    const int n_steps = 10000;         // Total steps
    const int output_freq = 100;       // Output every N steps
    const int thermo_freq = 1;
    const real tau_t = static_cast<real>(0.1f);
    
    // Box and lattice parameters
    const int n_cells = 5;             // 5x5x5 FCC lattice
    const real lattice_const = static_cast<real>(1.5f);  // Lattice constant
    const real box_size_1d = n_cells * lattice_const;
    Vector3 box_size(box_size_1d, box_size_1d, box_size_1d);
    
    // Create system
    const real cutoff = static_cast<real>(2.5f);         // LJ cutoff
    const real skin = static_cast<real>(0.3f);           // Neighbor list skin
    MDSystem system(box_size, dt, cutoff + skin);
    
    // Initialize particle positions on FCC lattice
    create_fcc_lattice(system, n_cells, lattice_const);
    std::cout << "Created " << system.n_particles() << " particles\n";
    std::cout << "  Ar: 80 particles\n";
    std::cout << "  Ne: 20 particles\n";
    
    // Initialize velocities
    system.initialize_velocities(temperature);
    std::vector<Vector3> vels = system.velocities();
    
    // Set up integrator and thermostat in the heap rather than 
    // stack to allow polymorphism throughout simulation steps
    auto integrator = std::make_unique<VelocityVerlet>();
    // system.set_integrator(std::move(integrator));
    
    auto thermostat = std::make_unique<BerendsenThermostat>(temperature, tau_t);
    // system.set_thermostat(std::move(thermostat));
    
    // Open output files
    std::ofstream energy_file("energy.dat");
    std::ofstream temp_file("temperature.dat");
    std::ofstream traj_file("trajectory.xyz");
    
    energy_file << "# Step  KineticE  PotentialE  TotalE\n";
    temp_file << "# Step  Temperature  Target_Temp\n";
    
    // Main simulation loop
    std::cout << "\nStarting simulation...\n";
    std::cout << "Step\tTemp\tKE\tPE\tTotal E\n";
    
    for (int step = 0; step <= n_steps; ++step) {
        // Compute forces
        system.compute_forces();
        
        // Apply thermostat periodically during equilibration
        if (step < n_steps/2 && step % thermo_freq == 0) {
            thermostat->apply(system);
        }
        
        // Integrate
        if (step > 0) {  // Skip first step integration
            integrator->step(system);
        }
        
        // Output
        if (step % output_freq == 0) {
            float ke = system.kinetic_energy();
            float pe = system.potential_energy();
            float total = ke + pe;
            float temp = system.temperature();
            
            // Console output
            std::cout << step << "\t" 
                     << std::fixed << std::setprecision(3)
                     << temp << "\t"
                     << ke << "\t"  
                     << pe << "\t"
                     << total << "\n";
            
            // File output
            energy_file << step << " " << ke << " " << pe << " " << total << "\n";
            temp_file << step << " " << temp << " " << temperature << "\n";
            
            // Write trajectory in XYZ format
            traj_file << system.n_particles() << "\n";
            traj_file << "Step " << step << " T=" << temp << "\n";
            
            const auto& positions = system.positions();
            const auto& types = system.types();
            
            for (size_t i = 0; i < system.n_particles(); ++i) {
                traj_file << types[i] << " " 
                         << positions[i].x() << " "
                         << positions[i].y() << " "  
                         << positions[i].z() << "\n";
            }
        }
    }
    
    // Final statistics
    std::cout << "\n=== Simulation Complete ===\n";
    std::cout << "Final temperature: " << system.temperature() << "\n";
    std::cout << "Final total energy: " << system.kinetic_energy() + system.potential_energy() << "\n";
    
    // Compute radial distribution function (optional)
    compute_rdf(system, "rdf.dat");
    
    energy_file.close();
    temp_file.close();
    traj_file.close();
    
    std::cout << "\nOutput files:\n";
    std::cout << "  energy.dat - Energy vs time\n";
    std::cout << "  temperature.dat - Temperature vs time\n";
    std::cout << "  trajectory.xyz - Trajectory (viewable in VMD/OVITO)\n";
    std::cout << "  rdf.dat - Radial distribution function\n";
    
    return 0;
}

// Optional: Compute radial distribution function
