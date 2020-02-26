#ifndef ODEINT_EIGEN_RK4_HPP
#define ODEINT_EIGEN_RK4_HPP
#include "../rk4.hpp"

namespace odeint
{

  using Vec_t = Eigen::VectorXd;
  using Mat_t = Eigen::MatrixXd;

  template<>
  RK4<Vec_t, Real>::RK4() {};

  template<>
  RK4<Vec_t, Real>::RK4(const Vec_t& initial_condition)
    : _x{initial_condition}, _x0{initial_condition} {};

  template<>
  RK4<Vec_t, Real>::RK4(const std::function<Vec_t(Real, const Vec_t&)>& f)
    : _eigen_rhs{f} {};

  template<>
  RK4<Vec_t, Real>::RK4(const std::function<PetscErrorCode(Real, const Vec_t&, Vec_t&)>& f) //wrong function type, so throw an exception
  {
    throw "Passed function with the signature required for RK4<Vec, Real> to RK4<Eigen::VectorXd, Real> ctor";
  }


  template<>
  RK4<Vec_t, Real>::RK4(const std::function<Vec_t(Real, const Vec_t&)>& f,
	const Vec_t& initial_condition, Real timestep)
    : _x{initial_condition}, _x0{initial_condition},
      _dt{timestep}, _eigen_rhs{f} {};

  template<>
  RK4<Vec_t, Real>::RK4(const std::function<PetscErrorCode(Real, const Vec_t&, Vec_t&)>& f,
	const Vec_t& initial_condition, Real timestep)
  {
    throw "Passed function with the signature required for RK4<Petsc::Vec, Real> to RK4<Vec_t, Real> ctor";
  }

  template<>
  RK4<Vec_t, Real>::RK4(const std::function<Vec_t(Real, const Vec_t&)>& f,
				  const Vec_t& initial_condition, Real dt,
				  Real time_window_length)
    : _x{initial_condition}, _x0{initial_condition},
      _dt{dt}, _tfinal{time_window_length}, _eigen_rhs{f} {};
  

  template<>
  RK4<Vec_t, Real>::RK4(const std::function<PetscErrorCode(Real, const Vec_t&, Vec_t&)>& f,
				  const Vec_t& initial_condition, Real timestep, Real time_window_length)
  {
    throw "Passed function with the signature required for RK4<Petsc::Vec, Real> to RK4<Eigen::VectorXd, Real> ctor";
  }

    template<>
  RK4<Vec_t, Real>::RK4(const std::function<Vec_t(Real, const Vec_t&)>& f,
				  const Vec_t& initial_condition, Real dt,
				  Size num_timesteps)
    : _x{initial_condition}, _x0{initial_condition},
      _dt{dt}, _tfinal{Real(num_timesteps) * dt},
      _num_timesteps{num_timesteps}, _eigen_rhs{f} {};

  template<>
  RK4<Vec_t, Real>::RK4(const std::function<PetscErrorCode(Real, const Vec_t&, Vec_t&)>& f,
				  const Vec_t& initial_condition, Real dt, Size num_timesteps)
  {
    throw "Passed function with the signature required for RK4<Petsc::Vec, Real> to RK4<Eigen::VectorXd, Real> ctor";
  }

  template<>
  void RK4<Vec_t, Real>::set_initial_condition(const Vec_t& x0)
  {
    _x0 = x0;
  }

  template<>
  void RK4<Vec_t, Real>::set_rhs(const std::function<Vec_t(Real, const Vec_t&)>& f)
  {
    _eigen_rhs = f;
  }

  template<>
  void RK4<Vec_t, Real>::set_rhs(const std::function<PetscErrorCode(Real, const Vec_t&, Vec_t&)>& f)
  {
    throw "Passed function with the signature required for RK4<Petsc::Vec, Real> to RK4<Eigen::VectorXd, Real>::set_rhs";
  }

  template<>
  RK4<Vec_t, Real>::~RK4() {};

  template<>
  void RK4<Vec_t, Real>::compute_k_values()
  {
    _k[0] = _dt * _eigen_rhs(_t, _x);
    _k[1] = _dt * _eigen_rhs(_t + 0.5*_dt, _x + 0.5 * _k[0]);
    _k[2] = _dt * _eigen_rhs(_t + 0.5*_dt, _x + 0.5 * _k[1]);
    _k[3] = _dt * _eigen_rhs(_t + _dt, _x + _k[2]);
  }

  template<>
  void RK4<Vec_t, Real>::advance_timestep()
  {
    if(_save_trajectory){
      _saved_trajectory.push_back({_t,_x});
    }
    compute_k_values();
    _x += _k[0] / 6.0 + _k[1] / 3.0 + _k[2] / 3.0 + _k[3] / 6.0;
    _t += _dt;
    ++_num_timesteps;
  }

  template<>
  Vec_t RK4<Vec_t, Real>::integrate()
  {
    _t = 0.0;
    _x = _x0;
    while(_t + _dt <= _tfinal){
      advance_timestep();
    }
    return _x;
  }

  template<>
  Vec_t RK4<Vec_t, Real>::integrate(Real time_length)
  {
    _tfinal = time_length;
    return integrate();
  }

  template<>
  Vec_t RK4<Vec_t, Real>::integrate(const Vec_t& x0)
  {
    _x0 = x0;
    return integrate();
  }

  
  template<>
  Vec_t RK4<Vec_t, Real>::integrate(const Vec_t& x0, Real time_length, bool save_trajectory)
  {
    _tfinal = time_length;
    _x0 = x0;
    _save_trajectory = save_trajectory;
    return integrate();
  }

  /*
  template<>
  Vec_t RK4<Vec_t, Real>::integrate(const Vec_t& x0, Size num_timesteps, bool save_trajectory)
  {
    set_num_timesteps(num_timesteps);
    _x0 = x0;
    _save_trajectory = save_trajectory;
    return integrate();
  }
  */

}//namespace odeint
#endif //ODEINT_EIGEN_RK4_HPP
