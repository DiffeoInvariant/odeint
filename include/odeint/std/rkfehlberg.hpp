#ifndef ODEINT_TEMPLATE_RKFEHLBERG_HPP
#define ODEINT_TEMPLATE_RKFEHLBERG_HPP
#include <cstddef>
#include <array>
#include <vector>
#include <utility>
#include <functional>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <string>
#include <omp.h>

namespace odeint
{

  static constexpr std::array<double, 5>   rkf_h_coeffs{0.25, 3.0/8.0,
						  12.0/13.0, 1.0, 0.5};
  static constexpr std::array<
      std::array<double, 5>, 5
    >                                      rkf_k_ccs{{
				 {0.25, 0, 0, 0, 0},
				 {(3.0/32.0), (9.0/32.0), 0 ,0, 0},
				 {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0, 0},
				 {439.0/216.0, -8.0, 3680.0/513.0, 845.0/4101.0, 0},
				 {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0}
						      }};
  static constexpr std::array<double, 4>   rk4_coeffs{25.0/216.0, 1408.0/2565.0,
						    2197.0/4101.0, -0.2};
  
  static constexpr std::array<double, 5>   rk5_coeffs{16.0/135.0, 6656.0/12825.0,
						  28561.0/56430.0, -9.0/50.0,
						  2.0/55.0};

  template<typename T, std::size_t N>
  static constexpr std::array<T, N> add(const std::array<T,N>& lhs, const std::array<T,N>& rhs) noexcept
  {
    auto res = lhs;
    for(auto i = 0; i < N; ++i){
      res[i] += rhs[i];
    }
    return res;
  }

  template<typename T, std::size_t N>
  static constexpr void inplace_add(std::array<T,N>& add_to, const std::array<T,N>& rhs) noexcept
  {
    for(auto i = 0; i < N; ++i){
      add_to[i] += rhs[i];
    }
  }

  template<typename T, std::size_t N, std::size_t M>
  static constexpr std::array<T,N> add_multiple(const std::array<std::array<T,N>, M>& arrays) noexcept
  {
    auto result = arrays[0];
    for(auto i = 1; i < M; ++i){
      inplace_add(result, arrays[i]);
    }
    return result;
  }
      

  template<typename T, std::size_t N>
  static constexpr std::array<T, N> scale(const std::array<T,N>& arr, T scale_by) noexcept
  {
    auto res = arr;
    for(auto i = 0; i < N; ++i){
      res[i] *= scale_by;
    }
    return res;
  }

  template<typename T, std::size_t N>
  static constexpr void inplace_scale(std::array<T,N>& arr, T scale_by) noexcept
  {
    for(auto i = 0; i < N; ++i){
      arr[i] *= scale_by;
    }
  }

  

  template<typename T, std::size_t N, typename Func>
  static std::array<std::array<T,N>,6> rkf_kvals(const Func& f, const std::array<T,N>& x,
				   T t, T dt) 
  {
    std::array<std::array<T,N>,6> rkf_k;
    rkf_k[0] = scale(f(t, x), dt);
    constexpr auto kcc1 = rkf_k_ccs[0];
    rkf_k[1] = scale(f(t + rkf_h_coeffs[0] * dt, add(x, scale(rkf_k[0], kcc1[0]))),
		     dt);


    constexpr auto kcc2 = rkf_k_ccs[1];
    rkf_k[2] = scale(f(t + rkf_h_coeffs[1] * dt,
		       add(x,
			   add(scale(rkf_k[1], kcc2[1]), scale(rkf_k[0], kcc2[0]))
			   )
		       ),
		     dt);
    
    constexpr auto kcc3 = rkf_k_ccs[2];
    rkf_k[3] = scale(f(t + rkf_h_coeffs[2] * dt,
		       add(x,
			   add(scale(rkf_k[2], kcc3[2]),
			       add(scale(rkf_k[1], kcc3[1]), scale(rkf_k[0], kcc3[0]))
			       )
			   )
		       ),
		     dt);

    constexpr auto kcc4 = rkf_k_ccs[3];
    rkf_k[4] = scale(f(t + rkf_h_coeffs[3] * dt,
		       add(x,
			   add(scale(rkf_k[3], kcc4[3]),
			       add(scale(rkf_k[2], kcc4[2]),
				   add(scale(rkf_k[1], kcc4[1]), scale(rkf_k[0], kcc4[0]))
				   )
			       )
			   )
		       ),
		     dt);

    constexpr auto kcc5 = rkf_k_ccs[4];
    rkf_k[5] = scale(f(t + rkf_h_coeffs[4] * dt,
		       add(x,
			   add(scale(rkf_k[4], kcc5[4]),
			       add(scale(rkf_k[3], kcc5[3]),
				   add(scale(rkf_k[2], kcc5[2]),
				       add(scale(rkf_k[1], kcc5[1]), scale(rkf_k[0], kcc5[0]))
				       )
				   )
			       )
			   )
		       ),
		     dt);

    return rkf_k;
  }

  template<typename T, std::size_t N>
  T norm(const std::array<T,N>& arr, int Norm=2)
  {
    
    if(Norm == -1){
      //use infinity norm
      auto max = std::abs(arr[0]);
      for(auto i = 1; i < N; ++i){
	if(std::abs(arr[i]) > max){
	  max = std::abs(arr[i]);
	}
      }
      return max;
    }

    T res;
    for(auto i = 0; i < N; ++i){
      res += std::pow(arr[i], Norm);
    }
    return std::pow(res, 1.0 / (static_cast<double>(Norm)));
  }
    
    


  template<typename T, std::size_t N>
  std::pair<std::array<T,N>,  T> rkfehlberg_deltas_and_resid(const std::array<std::array<T,N>,6>& rkf_k, T dt, int Norm=-1)
  {
    auto rk4delta = scale(rkf_k[0], rk4_coeffs[0]);
    inplace_add(rk4delta, scale(rkf_k[2], rk4_coeffs[1]));
    inplace_add(rk4delta, scale(rkf_k[3], rk4_coeffs[2]));
    inplace_add(rk4delta, scale(rkf_k[4], rk4_coeffs[3]));

    auto rk5delta = scale(rkf_k[0], rk5_coeffs[0]);
    inplace_add(rk5delta, scale(rkf_k[2], rk5_coeffs[1]));
    inplace_add(rk5delta, scale(rkf_k[3], rk5_coeffs[2]));
    inplace_add(rk5delta, scale(rkf_k[4], rk5_coeffs[3]));
    inplace_add(rk5delta, scale(rkf_k[5], rk5_coeffs[4]));
    
    return {rk4delta, norm(add(rk5delta, scale(rk4delta, -1.0)), Norm)};
  }

  //returns pair of (new value of dt, is resid < tol? (if yes use RK4 approx))
  template<typename T>
  std::pair<T, bool>
  recompute_rkfehlberg_dt(T resid, T tol, T t, T dt, T dt_min, T dt_max, T tfinal)
  {
    if(resid < tol){
      return {dt, true};
    }

    auto delta = 0.84 * std::pow(tol / resid, 0.25);
    if(delta <= 0.1){
      dt *= 0.1;
    } else if(delta >= 4.0){
      dt *= 4.0;
    } else {
      dt *= delta;
    }
    dt = std::min(dt, dt_max);
    if(t + dt > tfinal){
      dt = tfinal - t;
    } else if(dt < dt_min){
      throw std::runtime_error("Error, dt < minimum dt. Either the rhs function is too poorly-behaved to use Runge-Kutta-Fehlberg, or your dt_min is set too high.\n");
    }

    return {dt, false};
  }
    
      
    
    

  template<typename T, std::size_t N, typename Func>
  extern std::vector<std::pair<std::array<T, N>, T>>
  rkfehlberg_integrate(const Func& f, const std::array<T,N>& x0,
		       T dt_min, T dt_max, T t0, T tf, T tol, int TolNorm=-1)
  {
    std::vector<std::pair<std::array<T, N>, T>> states_at_times;
    //estimate max number of timesteps we'll take 
    states_at_times.reserve(static_cast<std::size_t>(2 * dt_min / (tf - t0)));
    auto t = t0;
    auto x = x0;
    auto dt = dt_max;
    states_at_times.push_back({x,t});
    while(t < tf){
      //std::cout << " time " << t << " and dt = " << dt << ".\n";
      auto kvals = rkf_kvals(f, x, t, dt);
      try {
      auto [rk4delta, resid] = rkfehlberg_deltas_and_resid(kvals, dt, TolNorm);
      auto [new_dt, use_rk4] = recompute_rkfehlberg_dt(resid, tol, t, dt, dt_min, dt_max, tf);
      dt = new_dt;
      if(use_rk4){
        inplace_add(x, rk4delta);
	t += dt;
	states_at_times.push_back({x,t});
      }
      } catch(const std::exception& e){
	throw;
      }
    }

    return states_at_times;
  }


  template<typename T, std::size_t N>
  void write_csv(std::string filename,
		 const std::vector<std::pair<std::array<T, N>, T>>& trajectory)
  {
    std::ofstream csv(filename);
    csv << "t";
    for(auto i = 0; i < N; ++i){
      csv << ",x_" << std::to_string(i);
    }
    csv << "\n";
    for(const auto& st : trajectory){
      auto state = st.first;
      auto time = st.second;
      csv << time;
      for(auto i = 0; i < N; ++i){
	csv << "," << state[i];
      }
      csv << "\n";
    }
  }
      
    
}
#endif //ODEINT_TEMPLATE_RKFEHLBERG_HPP
