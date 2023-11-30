#ifndef Newton_h
#define Newton_h

#include "lapack_interface.h"
#include "nonlinfunc.h"

namespace ASC_ode
{

  void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<double> res(func->DimF());
    Matrix<double, ColMajor> fprime(func->DimF(), func->DimX());

    for (int i = 0; i < maxsteps; i++)
      {
        func->Evaluate(x, res);
        // cout << "|res| = " << L2Norm(res) << endl;
        func->EvaluateDeriv(x, fprime);

        pep::bla::LapackLU tempLU(std::move(fprime));
        fprime = std::move(tempLU).Inverse();
        // CalcInverse(fprime);

        x -= fprime*res;
        double err= L2Norm(res);
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
