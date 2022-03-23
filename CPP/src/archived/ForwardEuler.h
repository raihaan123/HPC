#ifndef FORWARD_EULER_H
#define FORWARD_EULER_H

#include <fstream>

class ForwardEuler {
    protected:
        // Methods here can be accessed/overriden in derived classes! Simply forms the template.
        void TimeIntegrate();
        void SetParameters();
        void SetInitialConditions();

        // Virtual method provides a template for the derived classes to implement
        virtual void solve() = 0;

        std::fstream output;

    public:
        // Constructor
        ForwardEuler();

        // Destructor
        ~ForwardEuler();
};



#endif