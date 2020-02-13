#ifndef ODEINT_RK4_HPP
#define ODEINT_RK4_HPP

#include <omp.h>
#include "base.hpp"
#include <functional>
#include <optional>
#include <vector>
#include <array>

namespace odeint
{

  template<typename V=Eigen::VectorXd, typename T=Real>
  struct RK4
  {

    using StateType = V;
    
    RK4();

    RK4(const V& initial_condition);

    //takes the rhs function f, where x' = f(t, x). MUST have signature 
    // V(T, const V&) if V is an Eigen::VectorXd or similar, and
    // MUST have signature PetscErrorCode(T, const V&, V&) if V
    // is a Petsc Vec (the const V& is the current state y, the V& is an
    // output parameter)
    RK4(const std::function<V(T, const V&)>& f);

    //ctor for function taking Petsc Vecs
    RK4(const std::function<PetscErrorCode(T, const V&, V&)>& f);

    RK4(const std::function<V(T, const V&)>& f,
	const V& initial_condition);

    RK4(const std::function<PetscErrorCode(T, const V&, V&)>& f,
	const V& initial_condition);

    RK4(const std::function<V(T, const V&)>& f,
	const V& initial_condition, T timestep);

    RK4(const std::function<PetscErrorCode(T, const V&, V&)>& f,
	const V& initial_condition, T timestep);

    RK4(const std::function<V(T, const V&)>& f,
	const V& initial_condition, T dt,
	T time_window_length);

    RK4(const std::function<PetscErrorCode(T, const V&, V&)>& f,
	const V& initial_condition, T dt,
	T time_window_length);


    RK4(const std::function<V(T, const V&)>& f,
	const V& initial_condition, T dt,
        Size num_timesteps);

    RK4(const std::function<PetscErrorCode(T, const V&, V&)>& f,
	const V& initial_condition, T dt,
        Size num_timesteps);


    void set_initial_condition(const V& x0);

    V initial_condition()
    {
      return _x0;
    }

    V state()
    {
      return _x;
    }

    void set_final_time(T tfinal)
    {
      _tfinal = tfinal;
    }

    T final_time()
    {
      return _tfinal;
    }

    void set_num_timesteps(Size N)
    {
      _num_timesteps = N;
      _final_time = _dt * Real(N);
    }

    void set_dt(T new_dt)
    {
      _dt = new_dt;
      _num_timesteps = _tfinal / new_dt;
    }

    T dt()
    {
      return _dt;
    }

    void set_save_trajectory(bool save=true)
    {
      _save_trajectory = save;
    }

    std::optional<std::vector<V>> get_saved_trajectory()
    {
      if(_save_trajectory){
	return {_saved_trajectory};
      } else {
	return {};
      }
    }

    void get_petsc_options(PetscOptions opts_db=NULL)
    {
      if(has_petsc_option("-save_trajectory", std::nullopt, opts_db)){
	_save_trajectory = true;
      }
      auto [tfinal, has_tfinal] = get_petsc_option<Real>("-tf", std::nullopt, opts_db);
      if(has_tfinal){
	_tfinal = tfinal;
      }
      auto [dt, has_dt] = get_petsc_option<Real>("-dt", std::nullopt, opts_db);
      if(has_dt){
	_dt = dt;
      }
      auto [N, has_N] = get_petsc_option<Int>("-num_timesteps", std::nullopt, opts_db);
      if(has_N){
	_num_timesteps = static_cast<Size>(N);
      }

    }

    void set_rhs(const std::function<V(T, const V&)>& f);

    void set_rhs(const std::function<PetscErrorCode(T, const V&, V&)>& f);

    V integrate();

    V integrate(T time_length);

    V integrate(const V& x0);

    V integrate(const V& x0, T time_length, bool save_trajectory=false);

    V integrate(const V& x0, T dt, Size num_timesteps, bool save_trajectory=false);
    
    ~RK4();
    
  private:

    void compute_k_values();

    void advance_timestep();

    V                                  _x, _x0;
    T                                  _dt, _tfinal, _t = T(0.0);
    Size                               _num_timesteps=0;

    union
    {
    std::function<
      PetscErrorCode(T, const V&, V&)
    >                                  _petsc_rhs;

    std::function<
      V(T, const V&)
    >                                  _eigen_rhs;
    };
   
    bool                               _save_trajectory=false;
    std::vector<std::pair<T,V>>        _saved_trajectory;//pair of (time, state)

    std::array<V, 4>                   _k;

  };


}
#endif //ODEINT_RK4_HPP
