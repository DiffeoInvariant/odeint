#ifndef ODEINT_BASE_HPP
#define ODEINT_BASE_HPP
#ifndef ODEINT_NO_PETSC
#include <petscsys.h>
#include <petscvec.h>
#endif
#ifndef ODEINT_NO_EIGEN
#include <Eigen/Core>
#endif
#include <omp.h>
#include <string>
#include <cstddef>
#include <vector>
#include <optional>
#include <functional>
#include <utility>
namespace odeint
{




#ifndef ODEINT_NO_PETSC

  using Real = PetscReal;
  using Scalar = PetscScalar;
  using Size = std::size_t;
  using Int = PetscInt;
  
  PetscErrorCode initialize(int *argc, char ***argv,
		  const char file[]=NULL,
		  const char help[]=NULL)
  {
    auto ierr = PetscInitialize(argc, argv, file, help);CHKERRQ(ierr);
    return ierr;
  }


  bool has_petsc_option(std::string name,
			std::optional<std::string> prepend=std::nullopt,
			PetscOptions opts_db=NULL)
  {
    PetscBool has_opt;
    if(prepend){
      auto ierr = PetscOptionsHasName(opts_db, prepend->c_str(),
				      name.c_str(), &has_opt);CHKERRQ(ierr);
    } else {
      auto ierr = PetscOptionsHasName(opts_db, NULL,
				      name.c_str(), &has_opt);CHKERRQ(ierr);
    }
    return has_opt;
  }

  template<typename T>
  std::pair<T, bool> get_petsc_option(std::string name,
					    std::optional<std::string> prepend=std::nullopt,
					    PetscOptions opts_db=NULL);

  
  template<>
  std::pair<Real, bool> get_petsc_option<Real>(std::string name,
					       std::optional<std::string> prepend,
					       PetscOptions opts_db)
  {
    PetscBool has_opt;
    Real val=0;
    if(prepend){
      auto ierr = PetscOptionsGetReal(opts_db, prepend->c_str(),
				      name.c_str(), &val, &has_opt); 
    } else {
      auto ierr = PetscOptionsGetReal(opts_db, NULL,
				      name.c_str(), &val, &has_opt); 
    }
      
    return {val, has_opt};
  }

  template<>
  std::pair<Int, bool> get_petsc_option<Int>(std::string name,
					       std::optional<std::string> prepend,
					       PetscOptions opts_db)
  {
    PetscBool has_opt;
    Int val=0;
    if(prepend){
      auto ierr = PetscOptionsGetInt(opts_db, prepend->c_str(),
				     name.c_str(), &val, &has_opt); 
    } else {
      auto ierr = PetscOptionsGetInt(opts_db, NULL,
				     name.c_str(), &val, &has_opt); 
    }
      
    return {val, has_opt};
  }


  #define ODEINT_MAX_OPT_STR_LEN 100
  template<>
  std::pair<std::string, bool> get_petsc_option<std::string>(std::string name,
					                     std::optional<std::string> prepend,
							     PetscOptions opts_db)
  {
    PetscBool has_opt;
    char *val;

    if(prepend){
      auto ierr = PetscOptionsGetString(opts_db, prepend->c_str(),
				      name.c_str(), val, ODEINT_MAX_OPT_STR_LEN, &has_opt); 
    } else {
      auto ierr = PetscOptionsGetString(opts_db, NULL,
				      name.c_str(), val, ODEINT_MAX_OPT_STR_LEN, &has_opt); 
    }
    return {std::string{val}, has_opt};
  }


  #define ODEINT_MAX_OPT_ARR_LEN 10000
  template<>
  std::pair<std::vector<Int>, bool>
  get_petsc_option<std::vector<Int>>(std::string name,
				      std::optional<std::string> prepend,
				     PetscOptions opts_db)
  {
    PetscBool has_opt;
    Int size = ODEINT_MAX_OPT_ARR_LEN;
    Int *vals;
    if(prepend){
      auto ierr = PetscOptionsGetIntArray(opts_db, prepend->c_str(), name.c_str(),
					  vals, &size, &has_opt); 
    } else {
      auto ierr = PetscOptionsGetIntArray(opts_db, NULL, name.c_str(),
					  vals, &size, &has_opt); 
    }
      
    if(has_opt){
      return std::make_pair(std::vector(vals, vals+size), true);
    } else {
      return std::make_pair(std::vector<Int>(), false);
    }
  }

  template<>
  std::pair<std::vector<Real>, bool>
  get_petsc_option<std::vector<Real>>(std::string name,
				     std::optional<std::string> prepend,
				     PetscOptions opts_db)
  {
    PetscBool has_opt;
    Int size = ODEINT_MAX_OPT_ARR_LEN;
    Real *vals;
    if(prepend){
      auto ierr = PetscOptionsGetRealArray(opts_db, prepend->c_str(), name.c_str(),
					   vals, &size, &has_opt); 
    } else {
      auto ierr = PetscOptionsGetRealArray(opts_db, NULL, name.c_str(),
					   vals, &size, &has_opt); 
    }
    if(has_opt){
      auto retval = std::make_pair(std::vector(vals, vals+size), true);
      if(vals){
	delete[] vals;
	vals = NULL;
      }
      return retval;
    } else {
      return std::make_pair(std::vector<Real>(), false);
    }
  }

#else
  using Real = double;
  using Scalar = double;
  using Int = int;
  using Size = std::size_t;
  using PetscErrorCode = int;
#endif //ODEINT_NO_PETSC

  
  template<typename Func, typename... Args>
  std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>
  bind_rhs(const Func& f, Args&... args)
  {
    namespace ph = std::placeholders;
    return static_cast<
      std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>
      >(std::bind(f, ph::_1, ph::_2, args...));
  }
 




}//namespace odeint
#endif //ODEINT_BASE_HPP
