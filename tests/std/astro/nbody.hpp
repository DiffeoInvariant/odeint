#ifndef O_T_NBODY_H
#define O_T_NBODY_H
#include <odeint/std/rk4.hpp>
#include <odeint/std/rkfehlberg.hpp>
#include <cmath>

double norm_r_cubed(const std::array<double, 3>& r)
{
  double r3 = 0.0;
  for(const auto& rr : r){
    r3 += rr * rr;
  }
  return std::sqrt(r3) * std::sqrt(r3) * std::sqrt(r3);
}

class TwoBodyProblem
{
public:

  using Point = std::array<double, 3>;
  using State = std::array<double, 6>;
  TwoBodyProblem(const Point& init_r1, const Point& init_r2,
		 double _m1, double _m2, const Point& v0)
    : m1{_m1}, m2{_m2}, r0{init_r2}
  {
    for(auto i=0; i<3; ++i){
      r0[i] -= init_r1[i];
      x0[i] = r0[i];
    }
    for(auto j=3; j<6; ++j){
      x0[j] = v0[j-3];
    }
    
  };

  State initial_state() const noexcept
  {
    return x0;
  }

  auto operator()(double t, const State& x) const 
  {
    auto xprime = x;
    for(auto i=0; i<3; ++i){
      xprime[i] = x[3+i];
    }
    auto nr3 = norm_r_cubed({x[0],x[1],x[2]});
    for(auto i=3; i<6; ++i){
      xprime[i] = -G*(m1+m2) * x[i-3] / nr3;
    }
    return xprime;
  }
    


private:

  double m1, m2, G=1.0;
  Point  r0;
  State  x0;
};



#endif
