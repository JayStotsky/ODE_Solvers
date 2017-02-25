#ifndef _ODE_SOLVER_H_
#define _ODE_SOLVER_H_

//T is the return type i.e. double etc., S is the function type i.e. a function that returns a double - they must be compatible
template<class T, class S> class ODE_Solver
{
	public:
		ODE_Solver(){};
		~ODE_Solver(){};
		void GetOutput(Vector<double>&, Vector<T>& ) const;
		void SetInitialConditions(const T&);
		virtual void Solve(){};
		void SetFunction(std::shared_ptr<S>);  
		virtual void Initialize(){};
		bool CheckStatus() const;

	protected:
		bool T_start_set=0;
		bool T_end_set = 0;
		bool function_set=0;
		bool SolverDone = 0;
		
		double T_start=0;
		double T_end=0;
		Vector<T> Y;
		Vector<double> time;
		T Fty;
		std::shared_ptr<S> Function;
};

#include "ODE_SolverImplem.H"


#endif
