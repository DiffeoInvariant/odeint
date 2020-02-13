#include "odeint/std/rk4.hpp"
#include <iostream>
#include <fstream>
std::array<double, 2> forced_pendulum(double t, const std::array<double, 2>& x, double m, double l, double g, double beta, double A, double alpha)
{
  auto omega_prime = -beta * x[1] / m - g * std::sin(x[0])/l + A * std::cos(alpha * t) / (m * l);
  return {x[1], omega_prime};
}


int main()
{
  namespace ph = std::placeholders;
  

  double m = 0.1, l = 0.1, g = 9.8, beta = 0.25, A = 2.4, alpha = 1.59;
  auto rhs = std::bind(forced_pendulum, ph::_1, ph::_2, m, l, g, beta, A, alpha);

  std::array<double, 2> x0{1.0, -0.01};
  double t0 = 0.0, tf = 20.0, dt = 0.005;

  auto trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf);

  std::ofstream outfile("test_rk4_states.csv");

  for(const auto& st : trajectory){
    outfile << st.first[0] << ", " << st.first[1] << "\n";
  }

  return 0;
}
  
  



