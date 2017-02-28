#ifndef _NEWTONSOLVER_H_
#define _NEWTONSOLVER_H_

#include "LinearSolver.H"
#include "FunctionClasses.H"



template<class T> class NewtonSolver{
	public:
	NewtonSolver(){};
	~NewtonSolver(){};
	
	void SetJacobian(std::shared_ptr<JacobianFunction<T> >);
	void SetRHS(std::shared_ptr<RHSFunction<T> >);
	
	void SetTol(double);
	void SetMaxIter(int);
	void NewtonIteration();
	void Solve(T&, bool);
	
	private:
	LinearSolver<T,JacobianFunction<T> > L;
	std::shared_ptr<JacobianFunction<T> > dF;
	std::shared_ptr<RHSFunction<T> > F;
	
	T Y_old;
	T Y_new;
	T dY;
	T RHS;
	double tol=0.00001;
	double error=1.0;
	int max_iter = 100;
};


template<class T> class NewtonSolverTimeStep{
	public:
	NewtonSolverTimeStep(){};
	~NewtonSolverTimeStep(){};
	
	void SetTimeStepRoutine(std::shared_ptr<T> );	
	void SetTol(double);
	void SetMaxIter(int);
	void NewtonIteration();
	void Solve(T&, bool);
	
	private:
	T TimeStepRoutine;
	T Y_old;
	T Y_new;
	T dY;
	T RHS;
	double tol=0.00001;
	double error=1.0;
	int max_iter = 100;
};



#include "NewtonSolverImplem.H"


#endif
