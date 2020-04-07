CXX=clang++
CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -fopenmp

PETSC_DIR=/usr/local/petsc
PETSC_ARCH=arch-linux2-c-debug
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules
ODEINT_INCL=/home/diffeoinvariant/odeint/include/
INCLS=-I$(ODEINT_INCL) $(PETSC_CC_INCLUDES) `pkg-config --cflags eigen3`

TARGET=poincare_section

.PHONY: all clean $(TARGET) l63

all: $(TARGET) l63

l63: lorenz_gen.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(PETSC_WITH_EXTERNAL_LIB)

bbtraj: baseball_trajectory/batted_ball.cpp
	$(CXX) $(CXXFLAGS) $^ -o baseball_trajectory/distance $(INCLS) $(PETSC_WITH_EXTERNAL_LIB)

orbit: astro/nbody.cpp
	$(CXX) $(CXXFLAGS) $^ -o astro/orbit $(INCLS) $(PETSC_WITH_EXTERNAL_LIB)

$(TARGET): poincare_section_rk4.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(PETSC_WITH_EXTERNAL_LIB)

clean:
	@$(RM) *.o
