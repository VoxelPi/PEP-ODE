#include <nonlinfunc.h>
#include <ode.h>
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

int main() {
    {
    double tend = 4*M_PI;
    int steps = 100;
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

}
