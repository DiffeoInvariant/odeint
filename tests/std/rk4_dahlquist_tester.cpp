#include "odeint/std/rk4.hpp"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

using State = std::array<double, 1>;

State dahlquist(double t, const State& x, const double lambda)
{
  auto xprime = x;
  xprime[0] *= lambda;
  return xprime;
}

bool event_finder(double t, const State& x){
  return (x[0] <= 0.1 and x[0] >= 0.01);
}

void event_handler(double t, double& dt, const State& x,
		   odeint::RHSFunc_t<double, 1>& f)
{
  std::cout << "handling event: dt=" << dt << '\n';
  dt = 0.00001;
  std::cout << "handling event: dt=" << dt << '\n';
}

int main(){
  double t0 = 0.0, tf = 5.0, dt = 0.01, lambda = -0.8;
  State x0{1.0};
  namespace ph = std::placeholders;

  auto rhs = static_cast<odeint::RHSFunc_t<double, 1>>(std::bind(dahlquist, ph::_1, ph::_2, lambda));
  auto test_traj = odeint::rk4_integrate(rhs, x0, dt, t0, tf, {event_finder},
					 {event_handler});

  auto [xend, tend] = test_traj[test_traj.size() - 1];
  auto xe = xend[0];
  auto xendExact = std::exp(lambda * tend);

  std::cout << "Took " << test_traj.size() << " timesteps.\n";

  std::cout << "After 10 seconds integrating y' = -0.8  y, approximate solution is y = " << xe << ".\n Exact solution is " << xendExact << " and thus relative error is " << std::abs(xendExact - xe) / xendExact << ".\n";

  return 0;
}
  
