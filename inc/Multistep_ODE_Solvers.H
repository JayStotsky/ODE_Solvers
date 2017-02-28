#ifndef _MULTISTEP_ODE_SOLVERS_H_
#define _MULTISTEP_ODE_SOLVERS_H_

template<class T, class S> class ForwardEuler: public ODE_Solver<T,S>{
	public:
		ForwardEuler(){};
		~ForwardEuler(){};
		ForwardEuler(double tstart, double tend, int nsteps);
		ForwardEuler(double tstart, double tend, double dt);
		void Solve();
		void Solve(Vector<double>&, Vector<T>&);
		void Initialize(double,double,int);
		void Initialize(double,double,double);
		void StepForward();
		
	private:
		int Nsteps=0;
		double deltaT=0;
};

#include "Multistep_ODE_SolversImplem.H"

#endif 
