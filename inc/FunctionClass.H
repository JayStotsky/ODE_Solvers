#ifndef _FUNCTIONCLASS_H_
#define _FUNCTIONCLASS_H_

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"

template<class T> class ImplicitFixedStep;

template<class T> class ODELinearSolver
{
	public:
	void Solve(T&, T&, ImplicitFixedStep<T>*);
	
	void SetTol(double Tol)
	{ tol = Tol;};
	
	void SetMaxIter(int Max_iter)
	{max_iter = Max_iter;};
	
	private:
	double tol = 0.000001;
	double max_iter = 100;	
	T Ynew, Yold, Rnew, Rold, P, F;
	double alpha, beta;
	
};


template<class T> void ODELinearSolver<T>::Solve(T& Solution, T& RightHandSide, ImplicitFixedStep<T>* TimeStepper)
{
	int c = 0;
	Yold = Solution;
	Ynew = Solution;
	F = RightHandSide;
	TimeStepper->EvalJacobian(Ynew,F);
	Rold = RightHandSide-F;
	P = Rold;
	Rnew = Rold;
	double error=1.0;
	
	
	while (error>tol && c<max_iter)
	{
		TimeStepper->EvalJacobian(P,F);
		alpha = Dot(Rold,Rold)/Dot(P,F);
					
		Ynew = Yold+alpha*P;
		Rnew = Rold-alpha*F;
		 
		error = Norm(Rnew);
		beta = Dot(Rnew,Rnew)/Dot(Rold,Rold);
		
		P = Rnew+beta*P;
		
		Rold = Rnew;
		Yold = Ynew;
		c++;
	}
	if (c==max_iter && error>this->tol)
		printf("Warning: reached max number of iterations before error converged\n");
	Solution = Ynew;
	
};

template<class T> class ODENewton
{
	public:
	void Solve(T&, ImplicitFixedStep<T>*);
	
	void SetTol(double Tol)
	{tol = Tol;};
	
	void SetMaxIter(int Max_iter)
	{max_iter = Max_iter;};
		
	
	private:
		double tol = 0.000001;
		double max_iter = 100;
		T Yold, Ynew, Fy, Z;
		ODELinearSolver<T> LSolver;
};


template<class T> void ODENewton<T>::Solve(T& initial_guess, ImplicitFixedStep<T>* TimeStepper)
{
	int c = 0;
	double error=1.0;
	Yold = initial_guess;
	Ynew = initial_guess;
	Fy = initial_guess;
	Z = initial_guess;
	
	while(error>tol && c<max_iter)
	{
		TimeStepper->ConstructJacobian(Yold);
		TimeStepper->EvalRightHandSide(Yold,Fy);
		printf("%f\n", Fy);
		LSolver.Solve(Z,Fy,TimeStepper);
		Ynew = Yold + Z;
		error = Norm(Fy);
		Yold = Ynew;
		c++;
	}
	initial_guess = Ynew;
};

#endif
