#ifndef _ODE_SOLVER_H_
#define _ODE_SOLVER_H_

//T is the return type i.e. double etc., S is the function type i.e. a function that returns a double - they must be compatible

//Advantage of having class S as a template parameter: Generalizable to large range of possible functions
//Disadvantage of having class S as a template parameter: Can't reuse solver for two different functions that are different classes

//could create base function class and then define new functions through polymorphism?
template<class T> class ODE_Solver
{
	public:
		ODE_Solver(){};
		~ODE_Solver(){};
		void GetOutput(double&, T& ) const;
		void SetInitialConditions(const T&);
		virtual void Solve(){}; //don't print results
		virtual void Solve(Vector<double>&, Vector<T>&){}; //print results at every time step
		virtual void Solve(Vector<double>&, Vector<T>&, Vector<int>&){}; //print results at certain time steps
		virtual void StepForward(){}; //individual timestep function
		void SetFunction(std::shared_ptr<std::function<void(double&,T&,T&)> >);  
		virtual void Initialize(){};
		bool CheckStatus() const;

	protected:
		bool T_start_set=0;
		bool T_end_set = 0;
		bool function_set=0;
		bool SolverDone = 0;
		
		double T_start=0;
		double T_end=0;
		T Y;
		double time;
		T Fty;
		std::shared_ptr<std::function<void(double&,T&,T&)> > function;
};

#include "ODE_SolverImplem.H"


#endif
