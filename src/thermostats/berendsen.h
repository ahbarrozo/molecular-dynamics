#include <cmath>
#include <iostream>
#include <string>

#include "config.h"
#include "md_system.h"
#include "thermostats/thermostat.h"

class Berendsen : public Thermostat {
    private:
        real tau_;  // Coupling time constant
    public:
        Berendsen(real temp, real tau);
        void apply(MDSystem& system) override;
        std::string name() const override { return "Berendsen"; }
};