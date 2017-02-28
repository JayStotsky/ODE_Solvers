#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include "Vector.H"

template<class T> void ODE_Solver<T>::GetOutput(double& Time, T& Solution) const 
{
	if (SolverDone)
	{
		Time = time;
		Solution = Y;
	}
};

template<class T> bool ODE_Solver<T>::CheckStatus() const
{
	if (!T_start_set)
	{
		printf("Warning: Initial time not set, exiting solver\n");
		return 0;
	}	
	if (!T_end_set)
	{
		printf("Warning: Stopping time not set, exiting solver\n");
		return 0;
	}	
	if (!function_set)
	{
		printf("Warning: Function, f(t,y) not set, exiting solver\n");
		return 0;
	}	
	else
		return 1;
};

template<class T> void ODE_Solver<T>::SetInitialConditions(const T& Y0)
{
	if (CheckStatus())
	{
		time = T_start;
		Y = Y0;
		Fty = Y0;
		SolverDone = 0; //reset solver in case used previously
	}
	
};

template<class T> void ODE_Solver<T>::SetFunction(std::shared_ptr<std::function<void(double&,T&,T&)> > a_function)
{
	function = a_function;
	function_set = 1;
	SolverDone = 0;
};
