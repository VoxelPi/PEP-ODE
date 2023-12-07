#include <nonlinfunc.h>
#include <ode.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace pep::ode;

class MassSpring : public NonlinearFunction {
    public:
    size_t DimX() const override {
        return 2;
    }

    size_t DimF() const override {
        return 2;
    }
    
    void Evaluate (VectorView<double> x, VectorView<double> f) const override {
        f(0) = x(1);
        f(1) = -x(0);
    }
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override {
        df = 0.0;
        df(0, 1) = 1;
        df(1, 0) = -1;
    }
};

class ResistorCapacitor : public NonlinearFunction {
    public:
    double R = 1; // Ohm
    double C = 1e-3; // Farad

    size_t DimX() const override {
        return 2;
    }

    size_t DimF() const override {
        return 2;
    }
    
    void Evaluate (VectorView<double> x, VectorView<double> f) const override {
        f(0) = 1;
        f(1) = (std::cos(10.0 * M_PI * x(0)) - x(1)) / (R * C);
    }
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override {
        df = 0.0;
        df(0, 0) = 0;
        df(1, 1) = -1/(R*C);
        df(1, 0) = -1/(R*C) * 10 * M_PI * sin(10*M_PI*x(0));
    }
};

int main() {
    // Mass-Spring, implicit euler.
    {
        double tend = 4*M_PI;
        int steps = 200;
        Vector<double> y { 1, 0 };
        auto rhs = make_shared<MassSpring>();

        ofstream output;
        output.open ("implicit_mass_spring.csv");
        output << "t,y[0],y[1]\n";
        SolveODE_IE(tend, steps, y, rhs, [&output](double t, VectorView<double> y) {
            output << "" << t << "," << y(0) << "," << y(1) << std::endl;
        });
        output.close();
    }

    // Mass-Spring, explicit euler.
    {
        double tend = 4*M_PI;
        int steps = 200;
        Vector<double> y { 1, 0 };
        auto rhs = make_shared<MassSpring>();

        ofstream output;
        output.open ("explicit_mass_spring.csv");
        output << "t,y[0],y[1]\n";
        SolveODE_EE(tend, steps, y, rhs, [&output](double t, VectorView<double> y) {
            output << "" << t << "," << y(0) << "," << y(1) << std::endl;
        });
        output.close();
    }

    // Mass-Spring, crank nicolson.
    {
        double tend = 4*M_PI;
        int steps = 200;
        Vector<double> y { 1, 0 };
        auto rhs = make_shared<MassSpring>();

        ofstream output;
        output.open ("crank_nicolson_mass_spring.csv");
        output << "t,y[0],y[1]\n";
        SolveODE_CN(tend, steps, y, rhs, [&output](double t, VectorView<double> y) {
            output << "" << t << "," << y(0) << "," << y(1) << std::endl;
        });
        output.close();
    }

    // Mass-Spring, runge kutta.
    {
        // Matrix<double, ColMajor> a {
        //     { 0.0, 0.0, 0.0, 0.0 }, 
        //     { 0.5, 0.0, 0.0, 0.0 }, 
        //     { 0.0, 0.5, 0.0, 0.0 }, 
        //     { 0.0, 0.0, 1.0, 0.0 }, 
        // };
        // a = a.Transpose();
        // Vector<double> b { 1/6.0, 1/3.0, 1/3.0, 1/6.0 };
        // Vector<double> c { 0.0, 0.5, 0.5, 1.0 }; 
        Matrix<double, ColMajor> a {
            { 0.25, 0.25 - sqrt(3)/6 },
            { 0.25 + sqrt(3)/6, 0.25 },
        };
        Vector<double> b { 0.5, 0.5 };
        Vector<double> c { 0.5 - sqrt(3)/6, 0.5 + sqrt(3)/6 };
        // Matrix<double, ColMajor> a {
        //     { 5.0/12.0, -1.0/12.0 },
        //     { 3.0/4.0, 1.0/4.0 },
        // };
        // Vector<double> b { 0.75, 0.25 };
        // Vector<double> c { 1.0/3.0, 1.0 };

        double tend = 4*M_PI;
        int steps = 20;
        Vector<double> y { 1, 0 };
        auto rhs = make_shared<MassSpring>();

        ofstream output;
        output.open ("rk_mass_spring.csv");
        output << "t,y[0],y[1]\n";
        SolveODE_RK(tend, steps, a, b, c, y, rhs, [&output](double t, VectorView<double> y) {
            output << "" << t << "," << y(0) << "," << y(1) << std::endl;
        });
        output.close();
    }

    // Resistor-Capacitor, implicit euler.
    {
        double tend = 4/10.0;
        int steps = 100;
        Vector<double> y { 0, 0 };
        auto rhs = make_shared<ResistorCapacitor>();

        ofstream output;
        output.open ("implicit_rc.csv");
        output << "t,y[0],y[1]\n";
        SolveODE_IE(tend, steps, y, rhs, [&output](double t, VectorView<double> y) {
            output << "" << t << "," << y(0) << "," << y(1) << std::endl;
        });
        output.close();
    }

    // Resistor-Capacitor, explicit euler.
    {
        double tend = 4/10.0;
        int steps = 1000;
        Vector<double> y { 0, 0 };
        auto rhs = make_shared<ResistorCapacitor>();

        ofstream output;
        output.open ("explicit_rc.csv");
        output << "t,y[0],y[1]\n";
        SolveODE_EE(tend, steps, y, rhs, [&output](double t, VectorView<double> y) {
            output << "" << t << "," << y(0) << "," << y(1) << std::endl;
        });
        output.close();
    }

    // Resistor-Capacitor, crank nicolson.
    {
        double tend = 4/10.0;
        int steps = 250;
        Vector<double> y { 0, 0 };
        auto rhs = make_shared<ResistorCapacitor>();

        ofstream output;
        output.open ("crank_nicolson_rc.csv");
        output << "t,y[0],y[1]\n";
        SolveODE_EE(tend, steps, y, rhs, [&output](double t, VectorView<double> y) {
            output << "" << t << "," << y(0) << "," << y(1) << std::endl;
        });
        output.close();
    }
}
