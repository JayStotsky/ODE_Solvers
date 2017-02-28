//#include "NewtonSolver.H"
//#include "LinearSolver.H"
#include "FunctionClasses.H"
#include "Vector.H"

template<class T> void ImplicitFixedStep<T>::StepForward()
{
	if (JacobianSet)
	{
		N.Solve(this->Y,this);
		this->time = this->time+this->deltaT;
	}
	else
		printf("Error, implicit solver not set\n");
};
