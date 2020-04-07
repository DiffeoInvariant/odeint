#include "nbody.hpp"


int main(int argc, char **argv)
{

  using Point = std::array<double, 3>;

  Point r1{0,0,0};
  Point r2{1.0,0,0};
  Point vr2{0.0,1.0,0};
  double m1 = 0.5, m2 = 0.5, dt=0.001, tf=70.0;
  

  TwoBodyProblem tbp(r1, r2, m1, m2, vr2);

  //auto rhs = static_cast<odeint::RHSFunc_t<double, 6>>(tbp);
  auto trajectory = odeint::rk4_integrate(tbp, tbp.initial_state(), dt, 0.0, tf);

  odeint::write_csv("orbit1.csv", trajectory);

  return 0;
}
