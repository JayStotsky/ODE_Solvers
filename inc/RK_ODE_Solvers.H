#ifndef _RK_ODE_SOLVERS_H_
#define _RK_ODE_SOLVERS_H_

template<class T, class S> class RK4: public ODE_Solver<T,S>{
	public:
		RK4(){};
		~RK4(){};
		RK4(double tstart, double tend, int nsteps);
		RK4(double tstart, double tend, double dt);
		void Solve();
		void Solve(Vector<double>&, Vector<T>&);
		void StepForward();
		void Initialize(double,double,int);
		void Initialize(double,double,double);
	
	private:
		int Nsteps=0;
		double deltaT=0;
		T k1, k2, k3, k4;
		T f1, f2, f3;
		double t1, t2, t3;
};

#include "RK_ODE_SolversImplem.H"

#endif 
