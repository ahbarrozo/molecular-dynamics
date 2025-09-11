#pragma once
#include <string>
#include "integrators/integrator.h"

class VelocityVerlet : public Integrator {
    public:
        void step(MDSystem& system) override;
        std::string name() const override { return "VelocityVerlet"; }
};