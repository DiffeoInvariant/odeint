#ifndef ODEINT_TEMPLATE_RK4_HPP
#define ODEINT_TEMPLATE_RK4_HPP

#include <vector>
#include <array>
#include <cstddef>
#include <cmath>
#include <utility>
#include <functional>

namespace odeint
{
  inline namespace std
  {
    
  
  template<typename T, ::std::size_t N, typename Func>
  extern ::std::vector<::std::pair<::std::array<T, N>, T>>
  rk4_integrate(const Func& f, const ::std::array<T,N>& x0,
	        T dt, T t0, T tf)
  {
    ::std::vector<::std::pair<::std::array<T, N>, T>> results;
    auto num_timesteps = static_cast<::std::size_t>(::std::ceil((tf - t0)/dt));
    results.reserve(num_timesteps);
    
    auto x = x0;
    auto t = t0;
    results.push_back(::std::make_pair(x, t));
    for(auto i = 0; i < num_timesteps; ++i){
      auto k1 = f(t, x);
      auto xk = x;
      for(auto i = 0; i < N; ++i){
        xk[i] += 0.5 * k1[i];
      }
      
      auto k2 = f(t + 0.5 * dt, xk);
      for(auto i = 0; i < N; ++i){
        xk[i] = x[i] + 0.5 * k2[i];
      }
      auto k3 = f(t + 0.5 * dt, xk);
      for(auto i = 0; i < N; ++i){
        xk[i] = x[i] + k3[i];
      }
      auto k4 = f(t + dt, xk);

      for(auto i = 0; i < N; ++i){
        x[i] += dt * (k1[i] + k2[i] + k3[i] + k4[i]) / 6.0;
      }
      t += dt;
      results.push_back(::std::make_pair(x,t));
    }

    return results;
  }
      

    
  }//inline namespace std   

}//namespace odeint
#endif //ODEINT_TEMPLATE_RK4_HPP
