#include "odeint/adjoint.hpp"

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

struct Lorenz63RHS
{
  double a, b, r;

  Vec operator()(double t, const Vec& x)
  {
    auto xprime = x;
    xprime[0] = a * (x[1] - x[0]);
    xprime[1] = r * x[0] - x[1] - x[0] * x[2];
    xprime[2] = x[0] * x[1] - b * x[2];
    return xprime;
  }

};

struct Lorenz63RHSJacobian
{
  double a, b, r;

  void operator()(double t, const Vec& x, Mat& Jac)
  {
    Jac(0,0) = -a;
    Jac(0,1) = a;

    Jac(1,0) = r - x[2];
    Jac(1,1) = -1;
    Jac(1,2) = -x[0];

    Jac(2,0) = x[1];
    Jac(2,1) = x[0];
    Jac(2,2) = -b;
  }

};



int main(){

  double a = 16.0, b = 4.0, r = 45.0, dt = 0.001;
  
  std::size_t nstep = 100;
  Lorenz63RHS rhs{a,b,r};
  Lorenz63RHSJacobian rhs_jac{a,b,r};
  std::vector<Vec> x0s;
  Vec x01(3), x02(3), x03(3);
  x02 << 10.0, -5.0, 2.0;
  x01 << 0.0, 1.0, 2.0;
  x03 << 0.0, -1.0, 2.0;
  x0s.push_back(x01);
  x0s.push_back(x02);
  x0s.push_back(x03);

  for(const auto& x0 : x0s){
    odeint::RKAdjoint integrator(dt, nstep,x0, rhs, rhs_jac);
    auto [xfinal, dxfinaldx0] = integrator.integrate_with_grad();
    std::cout << "\n\nfinal values: \n" << xfinal << "\n";
    integrator.print_derivative(true);
  }

  return 0;
}
