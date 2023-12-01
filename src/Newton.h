#ifndef Newton_h
#define Newton_h

#include <iostream>
#include "lapack_interface.h"
#include "nonlinfunc.h"

namespace pep::ode {

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

        std::cout << "NEWTON ITERATION " << i << std::endl;
        std::cout << "    x = " << x << std::endl;
        std::cout << "    &x = " << &x << std::endl;
        std::cout << "    f(x) = " << res << std::endl;
        std::cout << "    &f(x) = " << &res << std::endl;
        std::cout << "    df(x) = " << fprime << std::endl;

        fprime =  pep::bla::LapackLU(fprime).Inverse();
        std::cout << "    df^-1 = " << fprime << std::endl;
        // CalcInverse(fprime);

        x -= fprime*res;
        std::cout << "    new x = " << x << std::endl;
        double err= L2Norm(res);
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
