#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>

#include "Vector.H"
#include "Matrix.H"
#include "FunctionClasses.H"

template<class T> class LinearSolver{
	public:
	LinearSolver(){};
	~LinearSolver(){};
	
	void SetMat(std::shared_ptr<std::function<void(T&,T&)> > M)
	{A = M;};
	void SetTol(double Tol)
	{tol = Tol;};
	void SetMaxIter(int Max_Iter)
	{max_iter = Max_Iter;};
	
	void SolveCG(T& Solution, T& RightHandSide)
	{
		int c = 0;
		x_old = RightHandSide;
		x_new = RightHandSide;
		A->eval(x_new, f);
		r_old = RightHandSide-f;
		p = r_old;
		while (error>tol && c<max_iter)
		{
			A->eval(p,f);
			alpha = Dot(r_old,r_old)/Dot(p,f);
						
			x_new = x_old+alpha*p;
			r_new = r_old-alpha*f;
			 
			error = Norm(r_new);
			beta = Dot(r_new,r_new)/Dot(r_old,r_old);
			
			p = r_new+beta*p;
			
			r_old = r_new;
			x_old = x_new;
			c++;
		}
		if (c==max_iter && error>tol)
			printf("Warning: reached max number of iterations before error converged\n");
		Solution = x_new;
	};
	
	void SolveMG(T& Solution, T& RightHandSide)
	{};
	
	void SolvePCG(T& Solution, T& RightHandSide)
	{};
	
	void SolveBCGStab(T& Solution, T& RightHandSide)
	{};
	
	void SolveAMG(T& Solution, T& RightHandSide)
	{};
	
	void SolveFAC(T& Solution, T& RightHandSide)
	{};
	
	void SolveLU(T& Solution, T& RightHandSide)
	{};
	
	void SolveCholesky(T& Solution, T& RightHandSide)
	{};
	
	
	private:
	std::shared_ptr<S> A;
	T x_old;
	T x_new;
	T r_old;
	T r_new;
	T p;
	T f;
	double tol=0.00001;
	double error=1.0;
	int max_iter = 100;
	double alpha, beta;
};

#endif
