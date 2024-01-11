#include <nonlinfunc.h>
#include <ode.h>
#include <iostream>
#include <fstream>

using namespace pep::ode;

class RHS : public NonlinearFunction
{
  size_t DimX() const override { return 1; }
  size_t DimF() const override { return 1; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = -x(0);
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
  {
    df(0, 0) = -1.0;
  }
};


int main()
{
  double tend = 10*M_PI;
  int steps = 1000;
  Vector<double> x { 1, };
  Vector<double> dx { 0. };
  auto rhs = make_shared<RHS>();
  auto mass = make_shared<IdentityFunction>(1);
  ofstream output;
  output.open ("mass_pendulum_newmark.csv");
  output << "t,x[0],x[1],x[2]\n";

  SolveODE_Newmark(
    tend, steps, x, dx, rhs, mass,
    [&output](double t, VectorView<double> x) {
      output << "" << t << "," << x(0) << std::endl;
    } 
      );
  output.close();
}
