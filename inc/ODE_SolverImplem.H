#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include "Vector.H"

template<class T, class S> void ODE_Solver<T,S>::GetOutput(double& Time, T& Solution) const 
{
	if (SolverDone)
	{
		Time = time;
		Solution = Y;
	}
};

template<class T, class S> bool ODE_Solver<T,S>::CheckStatus() const
{
	if (!T_start_set)
		return 0;
	if (!T_end_set)
		return 0;
	if (!function_set)
		return 0;
	else
		return 1;
};

template<class T, class S> void ODE_Solver<T,S>::SetInitialConditions(const T& Y0)
{
	if (CheckStatus())
	{
		time = T_start;
		Y = Y0;
		Fty = Y0;
	}
	
};

template<class T, class S> void ODE_Solver<T,S>::SetFunction(std::shared_ptr<S> function)
{
	Function = function;
	function_set = 1;
};
