#include <string>

#include "config.h"
#include "md_system.h"


class Thermostat {
    public:
        virtual ~Thermostat() = default;
        virtual void apply(MDSystem& system) = 0;
        virtual std::string name() const = 0;

        double get_target_temperature() const { return target_temperature_; }
        void set_target_temperature(double temp) { target_temperature_ = temp; }
    protected:
        float target_temperature_;
        explicit Thermostat(float temp) : target_temperature_(temp) {}
};