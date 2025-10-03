#include <cmath>
#include <iostream>
#include <string>

#include "config.h"
#include "md_system.h"
#include "thermostats/thermostat.h"

class BerendsenThermostat : public Thermostat {
    private:
        real tau_;  // Coupling time constant
    public:
        BerendsenThermostat(real temp, real tau);
        void apply(MDSystem& system) override;
        std::string name() const override { return "Berendsen"; }
};