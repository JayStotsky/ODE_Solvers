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
		void Initialize(double,double,int);
	
	private:
		int Nsteps=0;
		double deltaT=0;
};

#include "RK_ODE_SolversImplem.H"

#endif 
