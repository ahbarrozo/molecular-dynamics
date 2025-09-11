#pragma once

#include <memory>
#include <vector>

#include "md_system.h"
#include "potentials/potential.h"

class ForceField {
    private:
        // Non-bonded interactions
        std::vector<std::unique_ptr<Potential>> non_bonded_potentials_;
        
        // Bonded interactions (with connectivity)
        struct BondedInteraction {
            std::unique_ptr<Potential> potential;
            std::vector<size_t> particles;  // Indices of involved particles
        };
        
        std::vector<BondedInteraction> bonds_;
        std::vector<BondedInteraction> angles_;
        std::vector<BondedInteraction> dihedrals_;
        std::vector<BondedInteraction> impropers_;
        
        // Special interactions
        std::vector<std::unique_ptr<Potential>> external_potentials_;
        
    public:
        void add_non_bonded(std::unique_ptr<Potential> pot);
        void add_bond(std::unique_ptr<Potential> pot, size_t i, size_t j);
        void add_angle(std::unique_ptr<Potential> pot, size_t i, size_t j, size_t k);
        void add_dihedral(std::unique_ptr<Potential> pot, size_t i, size_t j, size_t k, size_t l);
        
        void compute_forces(MDSystem& system);
};