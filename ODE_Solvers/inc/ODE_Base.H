#ifndef _ODE_BASE_H_
#define _ODE_BASE_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"
#include "NonlinearSolver.H"

//makes sense to just have simple struct and then to have later on structs that are specialized, but will contain an instance of ODE_Info
struct ODE_Info
{
	double time;
	double dt;
	int nsteps ;
	bool defined;
	double start_time;
	double end_time;
};

template<typename T>
class ODE_Subclass
{
	public:
		ODE_Subclass(double,double,double);
		ODE_Subclass(double,double,int);
		virtual void RHSFunction(Vector<T>&,Vector<T>&){};
		virtual void LHSFunction(Vector<T>&,Vector<T>&){};
		void SetLHS(Vector<T>& y)
		{
			Yk = y;
		};
			
		
		void SetFunction(std::shared_ptr<std::function<void(double&,T&,T&)> >);
		void SetJacobian(std::shared_ptr<std::function<void(double&,T&,T&,T&)> >);
		void SetInitialCondition(T& Y0);
	
		bool ExplicitMethod = 0; 
		ODE_Info Info;
		Vector<T> Y; //current solution - vector since for with dG methods, solution is polynomial - i.e. multiple coefficients needed
		Vector<T> Yp; // working solution 
		Vector<T> Yk; // vector to specify in Jacobian
		Vector<T> Fyt; //right hand side vector
		
	protected:
		std::shared_ptr<std::function<void(double&,T&,T&)> > Function;
		std::shared_ptr<std::function<void(double&,T&,T&,T&)> > Jacobian;
		
};


template<typename T> 
class ODE_Superclass
{
	public:
		ODE_Superclass(std::shared_ptr<ODE_Subclass<T> >);
		
		void Solve(Vector<int>&, Vector<Vector<T> >&);
		void SolveTS_Implicit();
		void SolveTS_Explicit();
		void Print(int);
	
	private:
		std::shared_ptr<ODE_Subclass<T> > TS;	
		NonlinearSolver<T> NLS;	
		bool ExplicitMethod = 0;
		Vector<int>* store_time;
		Vector<Vector<T> >* store_Y;
		
};




#include "ODE_BaseImplem.H"

#endif
