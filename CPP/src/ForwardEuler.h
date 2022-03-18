#ifndef FORWARD_EULER_H
#define FORWARD_EULER_H


class ForwardEuler {
    protected:
        // Methods here can be accessed/overriden in derived classes! Simply forms the template.
        void TimeIntegrate();

    public:
        // Constructor
        ForwardEuler();
        void solve();

        // Destructor
        ~ForwardEuler();
};



#endif