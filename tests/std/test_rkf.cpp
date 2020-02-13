#include "odeint/std/rkfehlberg.hpp"
#include <functional>


using State = std::array<double, 3>;

State Lorenz63RHS(double t, const State& x, double a, double r, double b)
{
  State xprime = x;
  xprime[0] = a * (x[1] - x[0]);
  xprime[1] = r * x[0] - x[1] - x[0] * x[2];
  xprime[2] = x[0] * x[1] - b * x[2];
  return xprime;
}

State Rossler(double t, const State& x, double a, double b, double c)
{
  State xprime = x;
  xprime[0] = -(x[1] + x[2]);
  xprime[1] = x[0] + a * x[1];
  xprime[2] = b + x[2] * (x[0] - c);
  return xprime;
}

int main(int argc, char** argv)
{
  namespace ph = std::placeholders;
  
  const double tol = 1.0e-3, a = 16.0, r = 45.0, b = 4.0;
  const double dt_min=1.0e-5, dt_max=0.5, t0 = 0.0, tf = 30.0;
  //const State x0{-13.0, -12.0, 52.0};
  const std::vector<State> x0s{{ {-13.0, -12.0, 52.0}, {13.0, 12.0, -52.0},
				    {1.0, 5.0, 10.0}, {-30.2, -2.3, 1.24},
				    {24.0, 50.0, 130.0}
			      }};
  const std::vector<double> rvals{0.1, 0.3, 0.5, 0.7, 1.0, 5.0, 10.0,
				  13.5, 13.6, 14.7, 13.8, 13.9, 14.0,
				  20.0, 23.0, 23.5, 24.0, 24.5, 25.0,
				  25.5, 26.0, 26.5, 27.0, 27.5, 28.0,
				  28.5, 29.0, 29.5, 30.0, 40.0, 45.0,
				  60.0, 70.0, 75.0, 80.0, 90.0};
  #pragma omp parallel for
  {
    for(int j = 0; j < x0s.size(); ++j){
      auto x0 = x0s[j];
      for(int i = 0; i < rvals.size(); ++i){
    
	//std::cout << "trajectory with r = " << rvals[i] << ".\n";
	const auto rhs = std::bind(Lorenz63RHS, ph::_1, ph::_2, a, rvals[i], b);

	auto trajectory = odeint::rkfehlberg_integrate(rhs, x0, dt_min, dt_max, t0, tf, tol);

	auto filename = std::string{"test_rkf_lorenz63_r_"}
	                + std::to_string(rvals[i])
			+ std::string("_x0Num_")
			+ std::to_string(j+1)
			+ std::string{".csv"};
    
        odeint::write_csv(filename, trajectory);
      }
    }
  }

  return 0;
}
