CXX=clang++
CXXFLAGS=-std=c++17 -O3 -march=native -mtune=native -fopenmp

PETSC_DIR=/usr/local/petsc
PETSC_ARCH=arch-linux2-c-debug
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules
ODEINT_INCL=/home/diffeoinvariant/odeint/include/
INCLS=-I$(ODEINT_INCL) $(PETSC_CC_INCLUDES) `pkg-config --cflags eigen3`

TARGET=test_rk4

.PHONY: all clean $(TARGET)

all: $(TARGET)

$(TARGET): test_rk4.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(PETSC_WITH_EXTERNAL_LIB)

clean:
	@$(RM) *.o