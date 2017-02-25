CFLAGS = -g -Wall -std=c++11 -I. -I./inc -I./src
SRC = 
INC = Vector.H VectorImplem.H UtilityFunctions.H ODE_Solver.H ODE_SolverImplem.H Multipstep_ODE_Solvers.H Multistep_ODE_SolversImplem.H RK_ODE_Solvers.H RK_ODE_SolversImplem.H 
OBJ = 
CXX = g++

all: Test_ODE_Solver.exe
objects:

Test_ODE_Solver.exe: src/Test_ODE_Solver.cpp
	$(CXX) $(CFLAGS) $< -o $@
