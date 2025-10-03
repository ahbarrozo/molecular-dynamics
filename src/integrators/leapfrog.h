#pragma once
#include "integrators/integrator.h"
#include "md_system.h"

class Leapfrog : public Integrator {
private:
    bool first_step_ = true;  // Flag to handle first step differently
    
public:
    void reset() { first_step_ = true; }
    void step(MDSystem& system) override;
    std::vector<Vector3> get_synchronized_velocities(const MDSystem& system) const;
    std::string name() const override { return "Leapfrog"; }
};