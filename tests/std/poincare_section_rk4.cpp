#include "odeint/std/rk4.hpp"
#include "odeint/std/rkfehlberg.hpp"
#include "odeint/std/array_ops.hpp"
#include <cmath>

using State = std::array<double, 2>;

State forced_pendulum(double t, const State& x,
		      double m, double l, double g,
		      double beta, double A, double alpha)
{
  auto omega_prime = -beta * x[1] / m - g * std::sin(x[0])/l + A * std::cos(alpha * t) / (m * l);
  return {x[1], omega_prime};
}


struct PointCrossingFinder
{
  double T, dt;

  std::vector<std::pair<State, double>> crossings;

  PointCrossingFinder(double T0, double h) : T{T0}, dt{h}
  {};

  bool operator()(double t, const State& x)
  {
    auto res = std::fmod(t, T);
    return (res < dt);
  }

  void reset()
  {
    crossings.clear();
  }

};

struct InterpolatingCrossingFinder
{
  double T, dt, tprev;
  State xprev;
  std::vector<std::pair<State, double>> crossings;

  InterpolatingCrossingFinder(double T0, double h) : T{T0}, dt{h}
  {};

  bool operator()(double t, const State& x)
  {
    auto res = std::fmod(t, T);
    if(res >= T - dt){
      xprev = x;
      tprev = t;
    }

    return (res < dt);
  }

  void reset()
  {
    crossings.clear();
  }
};
    

void saveCrossing(PointCrossingFinder& finder, double t, double& dt, const State& x,
		    odeint::RHSFunc_t<double, 2>& f)
{
    finder.crossings.push_back({x, t});
}

void saveCrossing(InterpolatingCrossingFinder& finder, double t, double& dt, const State& x,
		    odeint::RHSFunc_t<double, 2>& f)
{
  double deltat = t - finder.tprev;
  auto n = std::floor(t / finder.T);
  double sec_time = n * finder.T;
  double dt_l = sec_time - finder.tprev;
  double wx = dt_l / deltat;
  double wprev = 1.0 - wx;

  auto interpolated_x = odeint::add(odeint::scale(finder.xprev, wprev), odeint::scale(x, wx));
  finder.crossings.push_back({interpolated_x, t});
}
     
      
  

int main(){

  namespace ph = std::placeholders;
  double m = 0.1, l = 0.1, g=9.8, beta = 0.0, alpha = 0.0, A = 0.0, t0 = 0.0, tf = 2000.0, dt=1.0e-4;
  State x0{0.01, 0.0};

  auto T0 = 2 * M_PI * std::sqrt(l / g);

  auto rhs = static_cast<odeint::RHSFunc_t<double, 2>>(std::bind(forced_pendulum, ph::_1, ph::_2, m, l, g, beta, A, alpha));
  PointCrossingFinder section_finder(T0, dt);

  auto save_crossing = [&section_finder](double t, double& dt, const State& x,
					 odeint::RHSFunc_t<double, 2>& f)
  {
    saveCrossing(section_finder, t, dt, x, f);
  };
  
  auto trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf, {section_finder},
					  {save_crossing});
  
  odeint::write_csv("p2a_section.csv", section_finder.crossings);

  section_finder.reset();

  section_finder.T = 0.42;

  trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf, {section_finder},
				     {save_crossing});

  odeint::write_csv("p2b_section_T_042.csv", section_finder.crossings);

  section_finder.reset();
  A = 1.8;
  beta = 0.25;
  alpha = 1.65;
  section_finder.T = 2 * M_PI / alpha;
  tf = section_finder.T * 2000;
  
  rhs = static_cast<odeint::RHSFunc_t<double, 2>>(std::bind(forced_pendulum, ph::_1, ph::_2, m, l, g, beta, A, alpha));
  x0 = {-1.4, -0.75};

  trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf, {section_finder},
				     {save_crossing});
  odeint::write_csv("p2c_traj_A18_b025_a165.csv", trajectory);
  odeint::write_csv("p2c_section_A18_b025_a165.csv", section_finder.crossings);


  InterpolatingCrossingFinder isection_finder(T0, dt);

  auto save_crossing_i = [&isection_finder](double t, double& dt, const State& x,
					 odeint::RHSFunc_t<double, 2>& f)
  {
    saveCrossing(isection_finder, t, dt, x, f);
  };

  beta = 0.0, alpha = 0.0, A = 0.0, t0 = 0.0, tf = 2000.0, dt=1.0e-4;
  x0 = {0.01, 0.0};
  rhs = static_cast<odeint::RHSFunc_t<double, 2>>(std::bind(forced_pendulum, ph::_1, ph::_2, m, l, g, beta, A, alpha));
  
  trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf, {isection_finder},
				     {save_crossing_i});
  
  odeint::write_csv("p3a_section.csv", isection_finder.crossings);

  isection_finder.reset();

  isection_finder.T = 0.87;

  trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf, {isection_finder},
				     {save_crossing_i});

  odeint::write_csv("p3b_section_T_087.csv", isection_finder.crossings);

  isection_finder.reset();
  A = 1.8;
  beta = 0.25;
  alpha = 1.65;
  isection_finder.T = 2 * M_PI / alpha;
  tf = isection_finder.T * 2000;

  rhs = static_cast<odeint::RHSFunc_t<double, 2>>(std::bind(forced_pendulum, ph::_1, ph::_2, m, l, g, beta, A, alpha));
  x0 = {-1.4, -0.75};

  trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf, {isection_finder},
				     {save_crossing_i});
  odeint::write_csv("p3c_traj_A18_b025_a165.csv", trajectory);
  odeint::write_csv("p3c_section_A18_b025_a165.csv", isection_finder.crossings);

  return 0;

}
  
