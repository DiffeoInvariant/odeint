#ifndef ODEINT_ARR_OPS_HPP
#define ODEINT_ARR_OPS_HPP
#include <array>

namespace odeint
{
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
    static constexpr std::array<T,N> add(const T *lhs_start, const T *rhs_start) 
    {
      std::array<T,N> result;
      for(auto i = 0; i < N; ++i){
	result[i] = lhs_start[i] + rhs_start[i];
      }

      return result;
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
  static constexpr void inplace_scale(std::array<T,N>& arr, const T scale_by) noexcept
  {
    for(auto i = 0; i < N; ++i){
      arr[i] *= scale_by;
    }
  }

    template<typename T>
    static constexpr void inplace_scale(T *start, T *end, const T scale_by) noexcept
    {
      for(; start != end; ++start){
	*start *= scale_by;
      }
    }

}

#endif
