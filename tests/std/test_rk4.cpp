#include "odeint/std/rk4.hpp"
#include "odeint/std/rkfehlberg.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>

std::array<double, 2> forced_pendulum(double t, const std::array<double, 2>& x, double m, double l, double g, double beta, double A, double alpha)
{
  auto omega_prime = -beta * x[1] / m - g * std::sin(x[0])/l + A * std::cos(alpha * t) / (m * l);
  return {x[1], omega_prime};
}

using State = std::array<double, 3>;

State Lorenz63RHS(double t, const State& x, double a, double r, double b)
{
  State xprime = x;
  xprime[0] = a * (x[1] - x[0]);
  xprime[1] = r * x[0] - x[1] - x[0] * x[2];
  xprime[2] = x[0] * x[1] - b * x[2];
  return xprime;
}


int main()
{
  namespace ph = std::placeholders;
  

  std::vector<double> betas{0.25};
  std::vector<double> alphas{/*0.05, 0.1, 0.2, 0.3, 0.4, 0.5,*/ 0.6, 0.7, 0.85, 1.0, 1.2, 1.3, 1.4, 1.5, 1.59, 1.60, 1.65};
  std::vector<double> As{0.5, 0.8, 1.0, 1.5, 1.8, 2.5};
  std::vector<std::array<double, 2>> x0s{{
					  /*{3.0, 0.1},
					  {0.01, 0.0},
					  {2.9, 0.1},
					  {2.7, 0.1},
					  {2.5, 0.1},
					  {2.3, 0.1},
					  {2.0, 0.1},
					  {1.7, 0.1},
					  {1.5, 0.1},
					  {1.3, 0.1},
					  {1.0, 0.1},
					  {0.8, 0.1},
					  {0.5, 0.1},
					  {0.4, 0.1},
					  {0.3, 0.1},
					  {0.2, 0.1},
					  {0.1, 0.1},
					  {0.0, 0.1},
					  {3.0, -0.1},
					  {2.9, -0.1},
					  {2.7, -0.1},
					  {2.5, -0.1},
					  {2.3, -0.1},
					  {2.0, -0.1},
					  {1.7, -0.1},
					  {1.5, -0.1},
					  {1.3, -0.1},
					  {1.0, -0.1},
					  {0.8, -0.1},
					  {0.5, -0.1},
					  {0.4, -0.1},
					  {0.3, -0.1},
					  {0.2, -0.1},
					  {0.1, -0.1},
					  {0.0, -0.1},
					  {3.0, 0.0},
					  {0.01, 0.0},
					  {2.9, 0.0},
					  {2.7, 0.0},
					  {2.5, 0.0},
					  {2.3, 0.0},
					  {2.0, 0.0},
					  {1.7, 0.0},
					  {1.5, 0.0},
					  {1.3, 0.0},
					  {1.0, 0.0},
					  {0.8, 0.0},
					  {0.5, 0.0},
					  {0.4, 0.0},
					  {0.3, 0.0},
					  {0.2, 0.0},
					  {0.1, 0.0},
					  {0.0, 0.0},
					  {-3.0, 0.0},
					  {-0.01, 0.0},
					  {-2.9, 0.0},
					  {-2.7, 0.0},
					  {-2.5, 0.0},
					  {-2.3, 0.0},
					  {-2.0, 0.0},
					  {-1.7, 0.0},
					  {-1.5, 0.0},
					  {-1.3, 0.0},
					  {-1.0, 0.0},
					  {-0.8, 0.0},
					  {-0.5, 0.0},
					  {-0.4, 0.0},
					  {-0.3, 0.0},
					  {-0.2, 0.0},
					  {-0.1, 0.0},
					  {-0.0, 0.0}*/

				 }};

  

  for(int j = 0; j < 30; ++j){
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), 0.5});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), -0.5});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), 1.0});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), -1.0});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), 0.25});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), -0.25});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), 0.75});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), -0.75});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), 1.5});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), -1.5});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), 2.0});
    x0s.push_back({3.0 - 0.2 * static_cast<double>(j), -2.0});
  }

  for(int j = 0; j < 10; ++j){
    x0s.push_back({0.9 + 0.2 * static_cast<double>(j), 0.5});
    x0s.push_back({-(0.9 + 0.2 * static_cast<double>(j)), -0.5});
    x0s.push_back({0.9 + 0.2 * static_cast<double>(j), 1.0});
    x0s.push_back({-(0.9 + 0.2 * static_cast<double>(j)), -1.0});
    x0s.push_back({0.9 + 0.2 * static_cast<double>(j), 0.25});
    x0s.push_back({-(0.9 + 0.2 * static_cast<double>(j)), -0.25});
    x0s.push_back({0.9 + 0.2 * static_cast<double>(j), 0.75});
    x0s.push_back({-(0.9 + 0.2 * static_cast<double>(j)), -0.75});
    x0s.push_back({0.9 + 0.2 * static_cast<double>(j), 1.5});
    x0s.push_back({-(0.9 + 0.2 * static_cast<double>(j)), -1.5});
    x0s.push_back({0.9 + 0.2 * static_cast<double>(j), 2.0});
    x0s.push_back({-(0.9 + 0.2 * static_cast<double>(j)), -2.0});
  }


  for(auto& x : x0s){
    if(std::abs(x[0]) < 0.15 and std::abs(x[1]) < 0.15){
      if(x[0] < 0){
	x[0] -= 0.5;
      } else {
	x[0] += 0.5;
      }
    }
  }
  
  double m = 0.1, l = 0.1, g = 9.8;
  double t0 = 0.0, tf = 10.0, dt = 0.001;

  for(const auto& beta : betas){
    for(const auto& alpha : alphas){
      for(const auto& A : As){
	
	auto rhs = std::bind(forced_pendulum, ph::_1, ph::_2, m, l, g, beta, A, alpha);
	#pragma omp parallel for
        for(int i = 0; i < x0s.size(); ++i){
	  auto x0 = x0s[i];
  
	auto trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf);
	/*
	for(auto i = 1; i < x0s.size(); ++i){
	  x0 = x0s[i];
	  auto new_traj = odeint::rk4_integrate(rhs, x0, dt, t0, tf);
	  trajectory.insert(trajectory.end(), new_traj.begin(), new_traj.end());
	  }*/

	  //State x0{-13.0,-12.0,52.0};
  //double t0 = 0.0, tf = 10.0, dt = 0.00001;
  //const auto rhs = std::bind(Lorenz63RHS, ph::_1, ph::_2, 16.0, 45.0, 4.0);
  //auto trajectory = odeint::rk4_integrate(rhs, x0, dt, t0, tf);
	
        std::string filename{"files/test_rk4_states"};
	filename += (std::string{"_alpha_"} + std::to_string(alpha)
		     + std::string{"_beta_"} + std::to_string(beta)
		     + std::string{"_A_"} + std::to_string(A)
		     + std::string{"_theta0_"} + std::to_string(x0[0])
		     + std::string{"_omega0_"} + std::to_string(x0[1])
		     + std::string{".csv"});
	

	/*for(const auto& st : trajectory){
	  outfile << st.first[0] << ", " << st.first[1] << ", " << st.first[2] << "\n";
	  }*/
	odeint::write_csv(filename, trajectory);
	}
      }
    }
  }

  return 0;
}
  
  



