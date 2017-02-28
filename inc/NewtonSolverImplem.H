#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>

#include "Vector.H"
#include "Matrix.H"
#include "LinearSolver.H"
#include "FunctionClasses.H"

template<class T> void NewtonSolver<T>::SetJacobian(std::shared_ptr<JacobianFunction<T> > M)
{
	dF = M;
	L.SetMat(dF);
};

template<class T> void NewtonSolver<T>::SetRHS(std::shared_ptr<RHSFunction<T> > R)
{
	F = R;
};

template<class T> void NewtonSolver<T>::SetTol(double Tol)
{
	tol = Tol;
};

template<class T> void NewtonSolver<T>::SetMaxIter(int Max_Iter)
{
	max_iter = Max_Iter;
};

template<class T> void NewtonSolver<T>::NewtonIteration()
{
	dY = Y_new;
	F->eval(Y_new,RHS);
	dF->setY(Y_old);
	L.SolveCG(dY,RHS);
	Y_new = Y_old+dY;
	Y_old = Y_new;
};

template<class T> void NewtonSolver<T>::Solve(T& Y, bool IsLinear)
{
	if (!IsLinear)
	{
		int c = 0;
		Y_new = Y;
		Y_old = Y;
		while (error>tol && c<max_iter)
		{
			NewtonIteration();
			error = Norm(Y_new-Y_old);
			c++;
		}
		Y = Y_new;
	}
	else
	{
		Y_new = Y;
		Y_old = Y;
		NewtonIteration();
		Y = Y_new;
	}
};


template<class T> void NewtonSolverTimeStep<T>::SetTimeStepRoutine(std::shared_ptr<T> TStep)
{
	TimeStepRoutine = TStep;
};

template<class T> void NewtonSolverTimeStep<T>::SetTol(double Tol)
{
	tol = Tol;
};

template<class T> void NewtonSolverTimeStep<T>::SetMaxIter(int Max_Iter)
{
	max_iter = Max_Iter;
};

template<class T> void NewtonSolverTimeStep<T>::NewtonIteration()
{
	dY = Y_new;
	TimeStepRoutine->RightHandSide(Y_new,RHS);
	TimeStepRoutine->SetJacobian(Y_old);
	TimeStepRoutine->LinearSolve(dY,RHS);
	Y_new = Y_old+dY;
	Y_old = Y_new;
};

template<class T> void NewtonSolverTimeStep<T>::Solve(T& Y, bool IsLinear)
{
	if (!IsLinear)
	{
		int c = 0;
		Y_new = Y;
		Y_old = Y;
		while (error>tol && c<max_iter)
		{
			NewtonIteration();
			error = Norm(Y_new-Y_old);
			c++;
		}
		Y = Y_new;
	}
	else
	{
		Y_new = Y;
		Y_old = Y;
		NewtonIteration();
		Y = Y_new;
	}
};
