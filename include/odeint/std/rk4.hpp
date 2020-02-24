#ifndef ODEINT_TEMPLATE_RK4_HPP
#define ODEINT_TEMPLATE_RK4_HPP

#include <vector>
#include <array>
#include <cstddef>
#include <cmath>
#include <utility>
#include <functional>
#include <iostream>
#include "array_ops.hpp"

namespace odeint
{
  inline namespace stdarray
  {

    std::array<double, 4>  _kcf{1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

    
    
  
  template<typename T, ::std::size_t N, typename Func>
  extern ::std::vector<::std::pair<::std::array<T, N>, T>>
  rk4_integrate(const Func& f, const ::std::array<T,N>& x0,
	        T dt, T t0, T tf)
  {
    std::array<std::array<T, N>, 5> _k;
    ::std::vector<::std::pair<::std::array<T, N>, T>> results;
    auto num_timesteps = static_cast<::std::size_t>(::std::ceil((tf - t0)/dt));
    results.reserve(num_timesteps);
    
    auto x = x0;
    auto t = t0;
    results.push_back(::std::make_pair(x, t));
    for(auto i = 0; i < num_timesteps; ++i){
      _k[0] = scale(f(t, x), dt);
      
      _k[1] = scale(f(t + 0.5 * dt, add(x, scale(_k[0], 0.5))), dt);
      _k[2] = scale(f(t + 0.5 * dt, add(x, scale(_k[1], 0.5))), dt);
      _k[3] = scale(f(t + dt, add(x, _k[2])), dt);
      //update
      _k[4] = scale(_k[0], _kcf[0]);

      for(auto j = 1; j < 4; ++j){
	inplace_add(_k[4], scale(_k[j], _kcf[j]));
      }

      inplace_add(x, _k[4]);
      t += dt;
      results.push_back(::std::make_pair(x,t));
    }

    return results;
  }

    template<typename T, std::size_t N>
    using RHSFunc_t = std::function<std::array<T,N>(T, const std::array<T,N>&)>;

  template<typename T, ::std::size_t N>
  extern ::std::vector<::std::pair<::std::array<T, N>, T>>
    rk4_integrate(RHSFunc_t<T,N>& f, const ::std::array<T,N>& x0,
	        T dt, T t0, T tf,
		const std::vector<
		std::function<bool(T, const std::array<T,N>&)>
		>& event_finders,
		const std::vector<
		  std::function<void(T, T&, const std::array<T,N>&, RHSFunc_t<T,N>&)>
		>& event_handlers)
  {
    if(event_finders.size() != event_handlers.size()){
	throw "Error, must supply the same number of event finders and event handlers.";
    }
    
    auto n_event_finders = event_finders.size();
    std::array<std::array<T, N>, 5> _k;
    ::std::vector<::std::pair<::std::array<T, N>, T>> results;
    auto num_timesteps = static_cast<::std::size_t>(::std::ceil((tf - t0)/dt));
    results.reserve(num_timesteps);
    
    auto x = x0;
    auto t = t0;
    results.push_back(::std::make_pair(x, t));
    for(auto i = 0; i < num_timesteps; ++i){
      for(auto ie = 0; ie < n_event_finders; ++ie){
	if((event_finders[ie])(t, x)){
	  (event_handlers[ie])(t, dt, x, f);
	}
      }
      _k[0] = scale(f(t, x), dt);
      
      _k[1] = scale(f(t + 0.5 * dt, add(x, scale(_k[0], 0.5))), dt);
      _k[2] = scale(f(t + 0.5 * dt, add(x, scale(_k[1], 0.5))), dt);
      _k[3] = scale(f(t + dt, add(x, _k[2])), dt);
      //update
      _k[4] = scale(_k[0], _kcf[0]);

      for(auto j = 1; j < 4; ++j){
	inplace_add(_k[4], scale(_k[j], _kcf[j]));
      }

      inplace_add(x, _k[4]);
      t += dt;
      results.push_back(::std::make_pair(x,t));
    }

    return results;
  }
      

    
  }//inline namespace std   

}//namespace odeint
#endif //ODEINT_TEMPLATE_RK4_HPP
