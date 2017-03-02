#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"
#include "ODE_Basic_Methods.H"
#include "ODE_Quadrature.H"
#include "Discrete_Galerkin_TS.H"

void scalarfunction(double& t, double &in, double &out)
{
	out = -in;
};

void scalarJacobian(double& t, double &y, double &in, double & out)
{
	out = -in;
}

void vectorfunction(double& t, Vector<double>& in, Vector<double>& out)
{	
	out = -in;

}

void vectorJacobian(double& t, Vector<double>& y, Vector<double>& in, Vector<double>& out)
{
	out = -in;
}

int main(int argc, char* argv[])
{
	double t0 = 0.0;
	double t1 = 4.0;
	int nsteps = 20;
	double dt = (t1-t0)/(double)nsteps;
	double y0 = 1.0;
	Vector<double> Y0(3);
	Y0(0) = 1.; Y0(1) = 1.; Y0(2) = 1.;
	Vector<int> t(nsteps+1);
	Vector<Vector<double> > y(nsteps+1);
	Vector<Vector<Vector<double> > > Y(nsteps+1);
	
	for (int i = 0; i < nsteps+1; i++)
	{
		t(i) = i;
		y(i).SetDim(1);
		Y(i).SetDim(3);
		Y(i)(0).SetDim(1);
	}
	//scalar forward Euler method test
	
	ForwardEuler<double> FE1(t0,t1,nsteps);
	FE1.SetInitialCondition(y0);
	FE1.SetFunction(std::make_shared<std::function<void(double&,double&,double&)> >(scalarfunction));
	std::shared_ptr<ForwardEuler<double> > FEptr = std::make_shared<ForwardEuler<double> >(FE1);
	ODE_Superclass<double> S2(FEptr);
	//S2.Solve(t,y);
	
		
	BackwardEuler<double> BE1(t0,t1,nsteps);
	BE1.SetInitialCondition(y0);
	BE1.SetFunction(std::make_shared<std::function<void(double&,double&,double&)> >(scalarfunction));
	BE1.SetJacobian(std::make_shared<std::function<void(double&,double&,double&,double&)> >(scalarJacobian));
	std::shared_ptr<BackwardEuler<double> > BEptr = std::make_shared<BackwardEuler<double> >(BE1);
	ODE_Superclass<double> S3(BEptr);
	S3.Solve(t,y);
			
		
	TrapezoidRule<double> T1(t0,t1,dt);
	T1.SetInitialCondition(y0);
	T1.SetFunction(std::make_shared<std::function<void(double&,double&,double&)> >(scalarfunction));
	T1.SetJacobian(std::make_shared<std::function<void(double&,double&,double&,double&)> >(scalarJacobian));
	std::shared_ptr<TrapezoidRule<double> > Tptr = std::make_shared<TrapezoidRule<double> >(T1);
	ODE_Superclass<double> S4(Tptr);
	//S4.Solve(t,y);
	
	
	//Vector test problems
	
	ForwardEuler<Vector<double> > FE2(t0,t1,dt);
	FE2.SetInitialCondition(Y0);
	FE2.SetFunction(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&)> >(vectorfunction));
	std::shared_ptr<ForwardEuler<Vector<double> > > FEptr2 = std::make_shared<ForwardEuler<Vector<double> > >(FE2);
	ODE_Superclass<Vector<double> > S5(FEptr2);
	//S5.Solve(t,Y);
	
	
	
	BackwardEuler<Vector<double> > BE2(t0,t1,0.01);
	BE2.SetInitialCondition(Y0);
	BE2.SetFunction(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&)> >(vectorfunction));
	BE2.SetJacobian(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&,Vector<double>&)> >(vectorJacobian));
	std::shared_ptr<BackwardEuler<Vector<double> > > BEptr2 = std::make_shared<BackwardEuler<Vector<double> > >(BE2);
	ODE_Superclass<Vector<double> > S6(BEptr2);
	//S6.Solve(t,Y);
	
	TrapezoidRule<Vector<double> > T2(t0,t1,0.01);
	T2.SetInitialCondition(Y0);
	T2.SetFunction(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&)> >(vectorfunction));
	T2.SetJacobian(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&,Vector<double>&)> >(vectorJacobian));
	std::shared_ptr<TrapezoidRule<Vector<double> > > Tptr2 = std::make_shared<TrapezoidRule<Vector<double> > >(T2);
	ODE_Superclass<Vector<double> > S7(Tptr2);
	//S7.Solve(t,Y);
	
	
	RungeKutta4<Vector<double> > RK4(t0,t1,dt);
	RK4.SetInitialCondition(Y0);
	RK4.SetFunction(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&)> >(vectorfunction));
	RK4.SetJacobian(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&,Vector<double>&)> >(vectorJacobian));
	std::shared_ptr<RungeKutta4<Vector<double> > > RKptr = std::make_shared<RungeKutta4<Vector<double> > >(RK4);
	ODE_Superclass<Vector<double> > S8(RKptr);
	S8.Solve(t,Y);
	
	
	MidpointRuleQuad<Vector<double>, DiscreteGalerkin0<Vector<double> > > Q(0,1);
	DiscreteGalerkin0<Vector<double> > DG0(t0,t1,dt);
	DG0.SetQuadrature(std::make_shared<ODE_Quadrature_Base<Vector<double>, DiscreteGalerkin0<Vector<double> > > >(Q));
	DG0.SetInitialCondition(Y0);
	DG0.SetFunction(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&)> >(vectorfunction));
	DG0.SetJacobian(std::make_shared<std::function<void(double&,Vector<double>&,Vector<double>&,Vector<double>&)> >(vectorJacobian));
	std::shared_ptr<DiscreteGalerkin0<Vector<double> > > DG0ptr = std::make_shared<DiscreteGalerkin0<Vector<double> > >(DG0);
	ODE_Superclass<Vector<double> > S9(DG0ptr);
	S9.Solve(t,Y);
	for (int i = 0; i < nsteps+1; i++)
		printf("%d %f %f\n", t(i), Y(i)(0)(0)-exp(-t(i)*dt), y(i)(0)-exp(-t(i)*dt)); 
	
};
