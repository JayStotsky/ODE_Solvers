#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include "Vector.H"
#include "ODE_Solver.H"
#include "Multistep_ODE_Solvers.H"
#include "RK_ODE_Solvers.H"
#include "TestClass.H"

class function{
	public:
		function(){};
		void eval(double& t, double& y, double &fty)
			{ fty = -y;};
	private:
};

class Function{
	public:
		Function(){};
		void eval(double& t, Vector<double>& y, Vector<double>& fty)
		{
			fty(0) = y(1);
			fty(1) = -y(0);
		}
			
	private:
};

int main(int argc, char* argv[])
{
	
	//test basic Forward Euler and Runge-Kutta 4th order method on scalar problem
	Vector<double> y;
	Vector<double> t;
	function f;
	std::shared_ptr<function> f_ptr = std::make_shared<function>(f);
	
	int Num_steps = 500;
	double T_start = 0;
	double T_end = 20;
	
	ForwardEuler<double,function> SolverFE(T_start,T_end, Num_steps);
	
	double y0 = 1;
	
    SolverFE.SetFunction(f_ptr);
	SolverFE.SetInitialConditions(y0);
	SolverFE.Solve(t,y);
	
	//printf("%d\n", y.GetDim());
	//for (int i =0; i<y.GetDim(); i++)
	//	printf("%f %f\n", t(i), y(i));


	RK4<double,function> SolverRK4(T_start,T_end,Num_steps);
	y0=1;
	SolverRK4.SetFunction(f_ptr);
	SolverRK4.SetInitialConditions(y0);
	SolverRK4.Solve(t,y);
	
	//printf("%d\n", y.GetDim());
	//for (int i =0; i<y.GetDim(); i++)
	//	printf("%f %f\n", t(i), y(i));
		
	//test basic Forward Euler and Runge-Kutta 4th order method on vector problem	
	ForwardEuler<Vector<double>, Function> SolverFE_Vector(T_start, T_end, Num_steps);
	Vector<Vector<double> > Y;
	
	Function F;
	std::shared_ptr<Function> F_ptr = std::make_shared<Function>(F);
	Vector<double> Y0(2); Y0(0) = 1; Y0(1) = 2;
	
	SolverFE_Vector.SetFunction(F_ptr);
	SolverFE_Vector.SetInitialConditions(Y0);
	SolverFE_Vector.Solve(t,Y);
	
	for (int i = 0; i < t.GetDim(); i++)
		printf("%f  %f  %f\n", t(i), Y(i)(0), Y(i)(1));
	

}

