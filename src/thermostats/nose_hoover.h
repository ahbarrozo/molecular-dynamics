#pragma once
#include "config.h"
#include "thermostat.h"

class NoseHooverThermostat : public Thermostat {
private:
    real tau_;
    real Q_;
    real xi_;
    real xi_dot_;
    
    // Chain parameters
    bool use_chains_;
    std::vector<real> xi_chain_;
    std::vector<real> xi_dot_chain_;
    int chain_length_;
    
public:
    // Simple Nosé-Hoover
    NoseHooverThermostat(real temp, real tau);
    
    // Nosé-Hoover chains (better ergodicity)
    NoseHooverThermostat(real temp, real tau, int chain_length);
    
    void apply(MDSystem& system) override;
    std::string name() const override { return "Nosè-Hoover"; }
    
    // Getters for monitoring
    real xi() const { return xi_; }
    real tau() const { return tau_; }
};