#ifndef _ODE_IMPLICITBASE_H_
#define _ODE_IMPLICITBASE_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"
#include "ODE_Base.H"


template<typename T>
class ODE_ImplicitBase
{
	public:
		ODE_ImplicitBase(double);
		ODE_ImplicitBase(int);
		void SetBounds(double, double);
		void SetInitialConditions(T &Y0);
		
		void SetODEFunction(std::shared_ptr<std::function<void(double&,T&,T&)> >);
		void SetJacobianFunction(std::shared_ptr<std::function<void(double&,T&,T&,T&)> >);
		
		void SetRHSFunction(std::shared_ptr<std::function<void(std::shared_ptr<std::function<void(double&,T&,T&)> >,T&, T&, ODE_Info*)> >);
		void SetLHSFunction(std::shared_ptr<std::function<void(std::shared_ptr<std::function<void(double&,T&,T&,T&)> >,T&, T&, ODE_Info*)> >);
		
		
		void StepForward();
		void Solve();
		
	protected:
		double start_time;
		double end_time;
		double time;
		double dt;
		int nsteps;
		
		ODE_Info Info;
		T Y;
		T Fy;
		
		std::shared_ptr<std::function<void(double&,T&,T&)> > ODEFunction;
		std::shared_ptr<std::function<void(double&,T&,T&,T&)> > Jacobian;
		std::shared_ptr<std::function<void(std::shared_ptr<std::function<void(double&,T&,T&)>) > >,T&, T&, ODE_Info*)> > RHS;
		std::shared_ptr<std::function<void(std::shared_ptr<std::function<void(double&,T&,T&,T&)>) > >,T&, T&, ODE_Info*)> > LHS;
};

#include "ODE_ImplicitBaseImplem.H"

#endif
