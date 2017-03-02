#ifndef _NONLINEAR_SOLVER_H_
#define _NONLINEAR_SOLVER_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"
#include "LinearSolver.H"

template<typename T> class ODE_Subclass; //forward declaration

template<typename T> class NonlinearSolver
{
	public:
		void SetODE(std::shared_ptr<ODE_Subclass<T> > ts)
		{
			TS = ts;
			L.SetODE(ts);
		};
		
		void Solve(Vector<T>& z)
		{
			Yold = TS->Y;
			Ynew = TS->Y;
			Fy = z; //initialize
			int c = 0; 
			double error = 1.0+tol;
			
			TS->RHSFunction(Yold,Fy);
			error = norm();
			Ynew = Yold;		
			while (error>tol && c<max_iter)
			{
				TS->SetLHS(Yold);
				L.Solve(Fy);
				Ynew = Fy+Yold;
				Yold = Ynew;	
				TS->RHSFunction(Yold,Fy);
				error = norm();	
				c++;
			}
			z = Ynew;			
		}
		double norm()
		{
			double error=0;
			for (int i = 0; i < Fy.GetDim(); i++)
				error+=pow(Norm(Fy(i)),2);
				
			return sqrt(error);
		}
		
	private:
		std::shared_ptr<ODE_Subclass<T> > TS;
		LinearSolver<T> L;
		double tol=0.000001;
		int max_iter=100;
		Vector<T> Yold, Ynew, Fy;
};








#endif
