#ifndef ODE_h
#define ODE_h

#include <iostream>
#include <functional>
#include <exception>

#include "Newton.h"

namespace pep::ode {
    
    void SolveODE_EE(
        double t_end,
        int steps,
        VectorView<double> y,
        shared_ptr<NonlinearFunction> rhs,
        std::function<void(double, VectorView<double>)> callback = nullptr
    ) {
        double dt = t_end / steps;
        Vector<double> res(rhs->DimF());

        double t = 0;
        for (int i = 0; i < steps; ++i) {
            t += dt;
            rhs->Evaluate(y, res);
            y += dt * res;
            if (callback) {
                callback(t, y);
            }
        }
    }

    // implicit Euler method for dy/dt = rhs(y)
    void SolveODE_IE(
        double tend, 
        int steps,
        VectorView<double> y, 
        shared_ptr<NonlinearFunction> rhs,
        std::function<void(double,VectorView<double>)> callback = nullptr
    ) {
        double dt = tend/steps;
        auto yold = make_shared<ConstantFunction>(y);
        auto ynew = make_shared<IdentityFunction>(y.Size());
        auto equ = ynew - yold - dt * rhs;

        double t = 0;
        for (int i = 0; i < steps; i++) {
            NewtonSolver(equ, y);
            yold->Set(y);
            t += dt;
            if (callback) {
                callback(t, y);
            }
        }
    }

    // Crank Nicolson method for dy/dt = rhs(y)
    void SolveODE_CN(
        double tend, 
        int steps,
        VectorView<double> y, 
        shared_ptr<NonlinearFunction> rhs,
        std::function<void(double,VectorView<double>)> callback = nullptr
    ) {
        double dt = tend/steps;
        Vector<double> res_old(rhs->DimF());
        rhs->Evaluate(y, res_old);

        auto yold = make_shared<ConstantFunction>(y);
        auto fold = make_shared<ConstantFunction>(res_old);
        auto ynew = make_shared<IdentityFunction>(y.Size());
        auto equ = ynew - yold - (dt/2.0) * (fold + rhs);

        double t = 0;
        for (int i = 0; i < steps; i++) {
            NewtonSolver(equ, y);

            rhs->Evaluate(y, res_old);
            fold->Set(res_old);
            yold->Set(y);
            t += dt;

            if (callback) {
                callback(t, y);
            }
        }
    }

    // Newmark and generalized alpha:
    // https://miaodi.github.io/finite%20element%20method/newmark-generalized/
    
    // Newmark method for  mass*d^2x/dt^2 = rhs
    void SolveODE_Newmark(
        double tend, 
        int steps,
        VectorView<double> x, VectorView<double> dx,
        shared_ptr<NonlinearFunction> rhs,   
        shared_ptr<NonlinearFunction> mass,  
        std::function<void(double,VectorView<double>)> callback = nullptr
    ) {
        double dt = tend/steps;
        double gamma = 0.5;
        double beta = 0.25;

        Vector<double> a(x.Size());
        Vector<double> v(x.Size());

        auto xold = make_shared<ConstantFunction>(x);
        auto vold = make_shared<ConstantFunction>(dx);
        auto aold = make_shared<ConstantFunction>(x);
        rhs->Evaluate (xold->Get(), aold->Get());
        
        auto anew = make_shared<IdentityFunction>(a.Size());
        auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
        auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

        auto equ = Compose(mass, anew) - Compose(rhs, xnew);

        double t = 0;
        for (int i = 0; i < steps; i++) {
            NewtonSolver (equ, a);
            xnew -> Evaluate (a, x);
            vnew -> Evaluate (a, v);

            xold->Set(x);
            vold->Set(v);
            aold->Set(a);
            t += dt;
            if (callback) callback(t, x);
        }
        dx = v;
    }

    // Generalized alpha method for M d^2x/dt^2 = rhs
    void SolveODE_Alpha (
        double tend, 
        int steps, 
        double rhoinf,
        VectorView<double> x, VectorView<double> dx, VectorView<double> ddx,
        shared_ptr<NonlinearFunction> rhs,   
        shared_ptr<NonlinearFunction> mass,  
        std::function<void(double,VectorView<double>)> callback = nullptr
    ) {
        double dt = tend/steps;
        double alpham = (2*rhoinf-1)/(rhoinf+1);
        double alphaf = rhoinf/(rhoinf+1);
        double gamma = 0.5-alpham+alphaf;
        double beta = 0.25 * (1-alpham+alphaf)*(1-alpham+alphaf);

        Vector<double> a(x.Size());
        Vector<double> v(x.Size());

        auto xold = make_shared<ConstantFunction>(x);
        auto vold = make_shared<ConstantFunction>(dx);
        auto aold = make_shared<ConstantFunction>(ddx);
        // rhs->Evaluate (xold->Get(), aold->Get()); // solve with M ???
        
        auto anew = make_shared<IdentityFunction>(a.Size());
        auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
        auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

        // auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - Compose(rhs, (1-alphaf)*xnew+alphaf*xold);
        auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - (1-alphaf)*Compose(rhs,xnew) - alphaf*Compose(rhs, xold);

        double t = 0;
        a = ddx;

        for (int i = 0; i < steps; i++) {
            NewtonSolver (equ, a);
            xnew -> Evaluate (a, x);
            vnew -> Evaluate (a, v);

            xold->Set(x);
            vold->Set(v);
            aold->Set(a);
            t += dt;
            if (callback) callback(t, x);
        }
        dx = v;
        ddx = a;
    }

    // implicit Euler method for dy/dt = rhs(y)
    void SolveODE_RK(
        double tend, 
        int steps,
        Matrix<double, ColMajor> a, 
        Vector<double> b, 
        Vector<double> c,
        VectorView<double> y, 
        shared_ptr<NonlinearFunction> rhs,
        std::function<void(double,VectorView<double>)> callback = nullptr
    ) {
        double dt = tend/steps;
        int s = c.Size();
        int n = y.Size();

        /*
        auto multiple_rhs = make_shared<MultipleFunc>(rhs, s);
        Vector<> my(s*y.Size());
        Vector<> mf(s*y.Size());  
        auto myold = make_shared<ConstantFunction>(my);
        auto mynew = make_shared<IdentityFunction>(s*n);
        auto equ = mynew-myold - dt * Compose(make_shared<MatVecFunc>(a, n), multiple_rhs);
                    

        double t = 0;
        for (int i = 0; i < steps; i++)
        {
        cout << "step " << i << endl;
        for (int j = 0; j < s; j++)
        my.Range(j*n, (j+1)*n) = y;
        myold->Set(my);

        NewtonSolver (equ, my);

        multiple_rhs->Evaluate(my, mf);
        for (int j = 0; j < s; j++)
        y += dt * b(j) * mf.Range(j*n, (j+1)*n);

        t += dt;
        if (callback) callback(t, y);
        }
        */

        // steps = 2;
        auto multiple_rhs = make_shared<BlockFunc>(rhs, s);
        Vector<double> mk(s*y.Size());
        Vector<double> my(s*y.Size());    
        auto myold = make_shared<ConstantFunction>(my);
        auto knew = make_shared<IdentityFunction>(s*n);
        auto equ = knew - Compose(multiple_rhs, myold+dt*make_shared<BlockMatVecFunc>(a, n));

        double t = 0;
        for (int i = 0; i < steps; i++) {
            // cout << "step " << i << endl;
            for (int j = 0; j < s; j++) {
                my.Range(j*n, (j+1)*n) = y;
            }
            myold->Set(my);

            mk = 0.0;
            NewtonSolver(equ, mk);

            for (int j = 0; j < s; j++) {
                y += dt * b(j) * mk.Range(j*n, (j+1)*n);
            }

            t += dt;
            if (callback) {
                callback(t, y);
            }
        }
    }
}

#endif
