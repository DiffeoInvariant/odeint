static char help[] = "Integrates Lorenz 63 for 0.1 seconds and gets the adjoint.";

#include <petscts.h>
#include <petsctao.h>

struct _lorenz_63_User {

  PetscReal a, b, r;
  PetscReal next_output, tf;

  PetscInt  nsteps;
  Mat       Jac;
  Vec       X, lambda[3];
};

typedef struct _lorenz_63_User *User;

static PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec X, Vec F, void *ctx)
{
  PetscErrorCode    ierr;
  User              user = (User)ctx;
  PetscScalar       *f;
  const PetscScalar *x;

  PetscFunctionBeginUser;
  
  ierr = VecGetArrayRead(X, &x);CHKERRQ(ierr);
  ierr = VecGetArray(F, &f);CHKERRQ(ierr);

  f[0] = user->a * (x[1] - x[0]);
  f[1] = user->r * x[0] - x[1] - x[0] * x[2];
  f[2] = x[0] * x[1] - user->b * x[2];

  ierr = VecRestoreArray(F, &f);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(X, &x);CHKERRQ(ierr);
				      
  PetscFunctionReturn(0);
}

static PetscErrorCode RHSJacobian(TS ts, PetscReal t, Vec X, Mat Jac, Mat Pre, void *ctx)
{
  PetscErrorCode    ierr;
  User              user = (User)ctx;
  PetscInt          rows_and_cols[] = {0,1,2};
  PetscScalar       J[3][3];
  const PetscScalar *x;

  PetscFunctionBeginUser;
  
  ierr = VecGetArrayRead(X, &x);CHKERRQ(ierr);
  J[0][0] = -1.0 * user->a;
  J[0][1] = user->a;
  J[0][2] = 0.0;

  J[1][0] = user->r - x[2];
  J[1][1] = -1.0;
  J[1][2] = -x[0];

  J[2][0] = x[1];
  J[2][1] = x[0];
  J[2][2] = -1.0 * user->b;

  ierr = MatSetValues(Jac, 3, rows_and_cols, 3, rows_and_cols, &J[0][0], INSERT_VALUES);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if(Pre != Jac){
    ierr = MatAssemblyBegin(Pre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Pre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X, &x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
  
 
  
int main(int argc, char **argv)
{
  TS ts;
  PetscMPIInt rank, size;
  struct _lorenz_63_User user;
  PetscReal dt = 0.001;
  PetscScalar *x0;/*[] = {10.0, -5.0, 2.0};*/
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help);CHKERRQ(ierr);
  user.a = 16.0;
  user.b = 4.0;
  user.r = 45.0;
  user.tf = dt * 100;
  user.nsteps = 0;

  ierr = MatCreate(PETSC_COMM_WORLD, &user.Jac);
  ierr = MatSetSizes(user.Jac, PETSC_DECIDE, PETSC_DECIDE, 3, 3);CHKERRQ(ierr);
  ierr = MatSetFromOptions(user.Jac);CHKERRQ(ierr);
  ierr = MatSetUp(user.Jac);CHKERRQ(ierr);
  ierr = MatCreateVecs(user.Jac, &user.X, NULL);CHKERRQ(ierr);

  ierr = TSCreate(PETSC_COMM_WORLD, &ts);CHKERRQ(ierr);
  ierr = TSSetEquationType(ts, TS_EQ_ODE_EXPLICIT);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts, NULL, RHSFunction, &user);CHKERRQ(ierr);
  ierr = TSSetRHSJacobian(ts, user.Jac, user.Jac, RHSJacobian, &user);CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts, user.tf);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts, dt);CHKERRQ(ierr);
  ierr = TSSetType(ts, TSRK);CHKERRQ(ierr);
  ierr = TSRKSetType(ts, TSRK4);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);

  ierr = VecGetArray(user.X, &x0);CHKERRQ(ierr);
  x0[0] = 10.0;
  x0[1] = -5.0;
  x0[2] = 2.0;
  ierr = VecRestoreArray(user.X, &x0);CHKERRQ(ierr);

  ierr = TSSetSaveTrajectory(ts);CHKERRQ(ierr);

  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ierr = TSSolve(ts, user.X);CHKERRQ(ierr);

  ierr = TSGetSolveTime(ts, &user.tf);CHKERRQ(ierr);
  ierr = TSGetStepNumber(ts, &user.nsteps);CHKERRQ(ierr);

  ierr = MatCreateVecs(user.Jac, &user.lambda[0], NULL);
  ierr = VecGetArray(user.lambda[0], &x0);CHKERRQ(ierr);
  x0[0] = 1.0; x0[1] = 0.0; x0[2] = 0.0;
  ierr = VecRestoreArray(user.lambda[0], &x0);CHKERRQ(ierr);

  ierr = MatCreateVecs(user.Jac, &user.lambda[1], NULL);
  ierr = VecGetArray(user.lambda[1], &x0);CHKERRQ(ierr);
  x0[0] = 0.0; x0[1] = 1.0; x0[2] = 0.0;
  ierr = VecRestoreArray(user.lambda[1], &x0);CHKERRQ(ierr);

  ierr = MatCreateVecs(user.Jac, &user.lambda[2], NULL);
  ierr = VecGetArray(user.lambda[2], &x0);CHKERRQ(ierr);
  x0[0] = 0.0; x0[1] = 0.0; x0[2] = 1.0;
  ierr = VecRestoreArray(user.lambda[2], &x0);CHKERRQ(ierr);

  ierr = TSSetCostGradients(ts, 3, user.lambda, NULL);CHKERRQ(ierr);

  ierr = TSAdjointSolve(ts);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Derivative rows (adjoint variables):\n");CHKERRQ(ierr);
  for(auto i = 0; i < 3; ++i){
    ierr = VecView(user.lambda[i], PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");CHKERRQ(ierr);
  }

  ierr = MatDestroy(&user.Jac);CHKERRQ(ierr);
  ierr = VecDestroy(&user.X);CHKERRQ(ierr);
  ierr = VecDestroy(&user.lambda[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&user.lambda[1]);CHKERRQ(ierr);
  ierr = VecDestroy(&user.lambda[2]);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return ierr;
     
}
