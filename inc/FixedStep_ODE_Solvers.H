#ifndef _FIXEDSTEP_ODE_SOLVERS_H_
#define _FIXEDSTEP_ODE_SOLVERS_H_
#include "ODE_Solver.H"

template<class T> class FixedStep_ODE_Solver: public ODE_Solver<T>{
	public:
		FixedStep_ODE_Solver(){};
		~FixedStep_ODE_Solver(){};
		FixedStep_ODE_Solver(double tstart, double tend, int nsteps);
		FixedStep_ODE_Solver(double tstart, double tend, double dt);
		void Solve();
		void Solve(Vector<double>&, Vector<T>&);
		void Initialize(double,double,int);
		void Initialize(double,double,double);
		virtual void StepForward(){};
		
	protected:
		int Nsteps=0;
		double deltaT=0;
		bool ready = 0;
};

#include "FixedStep_ODE_SolversImplem.H"

#endif 
