#include <string>

#include "md_system.h"

class Integrator {
    public:
        virtual ~Integrator() = default;
        virtual void step(MDSystem& system) = 0;
        virtual void reset() = 0;
        virtual std::string name() const = 0;
};