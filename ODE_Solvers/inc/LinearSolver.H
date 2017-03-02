#ifndef _LINEAR_SOLVER_H_
#define _LINEAR_SOLVER_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"

template<typename T> class ODE_Subclass; //forward declaration

template<typename T> class LinearSolver
{
	public:
		void SetODE(std::shared_ptr<ODE_Subclass<T> > ts)
		{
			TS = ts;
		};
		
		void Solve(Vector<T>& Ynew)
		{
			//Ynew serves both as the initial right hand side, and then as the output - easy enough to change this using a two argument function if later needed
			int c = 0; 
			double error = 1.0+tol;
			
			AX = TS->Y; //initialize vector
			TS->LHSFunction(Ynew,AX);
			Rold = Ynew - AX;  
			P = Rold;
			Xold = Ynew;
			Xnew = Ynew;
			Rnew = Rold;
			
			
			while (error>tol && c<max_iter)
			{
				TS->LHSFunction(P,AX);
				alpha = dot(Rold,Rold)/dot(P,AX);
				
				for (int i = 0;  i < P.GetDim(); i++)
				{	
					Xnew(i) = Xold(i) + alpha * P(i);
					Rnew(i) = Rold(i) - alpha * AX(i);
				}
				error = norm();
				
				beta = dot(Rnew, Rnew)/dot(Rold,Rold);
				for (int i = 0; i < P.GetDim(); i++)
					P(i) = Rnew(i) + beta * P(i);
				c++;
				//linear solver iterations
			}
			Ynew = Xnew;
			
			
			//or direct solver
		}
		
		//need norm and dot product of potential "vector of vectors"
		double norm()
		{
			double error = 0;
			for (int i = 0; i<Rnew.GetDim(); i++)
				error+=pow(Norm(Rnew(i)),2);
			return sqrt(error);
		};
		double dot(Vector<T>& V1, Vector<T>& V2)
		{
			double dotprod = 0;
			for (int i = 0; i < V1.GetDim(); i++)
				dotprod+=Dot(V1(i),V2(i));
			return dotprod;
		};
		
	private:
		std::shared_ptr<ODE_Subclass<T> > TS;
		double tol=0.00001;
		int max_iter=100;
		Vector<T> Xold, Xnew, Rold, Rnew, P, AX;
		double alpha, beta;
};








#endif
