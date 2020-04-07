#ifndef ODEINT_ADJOINT_HPP
#define ODEINT_ADJOINT_HPP
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <functional>
#include <utility>
#include <vector>
#include <array>
#include <iostream>

namespace odeint
{

  struct RKAdjoint
  {
    using Vec_t = Eigen::VectorXd;
    using Mat_t = Eigen::MatrixXd;

    RKAdjoint(double h, std::size_t num_steps, const Vec_t& x_0,
	      const std::function<Vec_t(double, const Vec_t&)>& rhs,
	      const std::function<void(double, const Vec_t&, Mat_t&)> rhs_jacobian,
	      bool save_traj=false)
      : t{0.0}, t0{0.0}, tf{static_cast<double>(h * num_steps)}, dt{h}, nstep{num_steps}, x{x_0}, x0{x_0},
	Jac(Mat_t::Zero(x_0.size(), x_0.size())),
	delta(Mat_t::Identity(x_0.size(), x_0.size())),
	delta0(Mat_t::Identity(x_0.size(), x_0.size())),
	rhs_jac{rhs_jacobian},
	rhs_func{rhs}, N{x_0.size()}, save{save_traj}
    {};

    //derivative of final value wrt initial condition
    Mat_t getDerivative()
    {
      return delta;
    }

    //returns column sums of the derivative
    Vec_t getVariations()
    {
      auto ones = Vec_t::Ones(N);
      return ones.transpose() * delta;
    }
	      
    Vec_t integrate()
    {
      for(auto i = 0; i < nstep; ++i){
	advance_timestep();
      }
      return x;
    }

    Eigen::VectorXcd lyapunovExponents()
    {
      auto ev = delta.eigenvalues();
      for(auto i=0; i<ev.size(); ++i){
	ev[i] = std::log(std::abs(ev[i]))/t;
      }
      return ev;
    }

    std::pair<Vec_t, Mat_t> integrate_with_grad()
    {
      for(auto i = 0; i < nstep; ++i){
	advance_timestep();
      }
      return std::make_pair(x, delta);
    }

    void print_derivative(bool with_variations=false)
    {
      int io_precision = 5;
      Eigen::IOFormat clean(io_precision, 0, ",","\n","[","]");
      std::cout << "Derivative w.r.t. initial condition:\n" << delta.format(clean) << '\n';

      if(with_variations){
	auto var = getVariations();
	std::cout << "Variations after evolution:\n" << var.format(clean) << '\n';
      }
    }
			    
      

  private:

    inline void compute_k_values()
    {
      k[0] = dt * rhs_func(t, x);
      k[1] = dt * rhs_func(t + 0.5 * dt, x + 0.5 * k[0]);
      k[2] = dt * rhs_func(t + 0.5 * dt, x + 0.5 * k[1]);
      k[3] = dt * rhs_func(t + dt, x + k[2]);
    }

    inline void update_jacobian(double tt, const Vec_t& xx)
    {
      rhs_jac(tt, xx, Jac);
      //std::cout << "jacobian \n" << Jac << "\n eigenvalues \n" << Jac.eigenvalues() << "\n\n";
    }

    inline void compute_delta_k_values()
    {
      update_jacobian(t, x);
      dk[0] = dt * Jac * delta;
      update_jacobian(t + 0.5 * dt, x + 0.5 * k[0]);
      dk[1] = dt * Jac * (delta + 0.5 * dk[0]);
      update_jacobian(t + 0.5 * dt, x + 0.5 * k[1]);
      dk[2] = dt * Jac * (delta + 0.5 * dk[1]);
      update_jacobian(t + dt, x + k[2]);
      dk[3] = dt * Jac * (delta + dk[2]);
    }

    inline void advance_delta()
    {
      delta += dk[0] / 6.0 + dk[1] / 3.0 + dk[2] / 3.0 + dk[3] / 6.0;
    }

    inline void advance_x_t()
    {
      x += k[0] / 6.0 + k[1] / 3.0 + k[2] / 3.0 + k[3] / 6.0;
      t += dt;
    }

    inline void advance_timestep()
    {
      if(save){
	t_x_stored.push_back({t, x});
      }
      compute_k_values();
      compute_delta_k_values();
      advance_delta();
      advance_x_t();
    }
    
      
    std::array<Vec_t, 4> k;
    std::array<Mat_t, 4> dk;
    double             t, t0, tf, dt;
    std::size_t        nstep;
    Vec_t              x, x0;
    Mat_t              Jac, delta, delta0;
    std::function<void(double, const Vec_t&, Mat_t&)> rhs_jac;
    std::function<Vec_t(double, const Vec_t&)> rhs_func;
    long int                N;
    bool               save;
    std::vector<std::pair<double, Vec_t>> t_x_stored;
  };







}

#endif
