#ifndef _ODE_BASIC_METHODS_H_
#define _ODE_BASIC_METHODS_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"
#include "ODE_Base.H"

template<typename T>
class ForwardEuler: public ODE_Subclass<T>
{
	public:
		ForwardEuler(double t0,double t1,double dt):ODE_Subclass<T>(t0,t1,dt)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Fyt.SetDim(1);
			this->ExplicitMethod = 1;
		};
		ForwardEuler(double t0,double t1,int nsteps):ODE_Subclass<T>(t0,t1,nsteps)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->ExplicitMethod = 1;
		};
		void RHSFunction(Vector<T>& y,Vector<T>& yp)
		{
			this->Fyt = y;
			(*(this->Function))( (this->Info.time), y(0), this->Fyt(0));
			yp(0) = y(0) + (this->Info.dt)*(this->Fyt(0));
		};
};


template<typename T>
class BackwardEuler: public ODE_Subclass<T>
{
	public:
		BackwardEuler(double t0,double t1,double dt):ODE_Subclass<T>(t0,t1,dt)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Fyt.SetDim(1);
			this->ExplicitMethod = 0;
		};
		BackwardEuler(double t0,double t1,int nsteps):ODE_Subclass<T>(t0,t1,nsteps)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->ExplicitMethod = 0;
		};
		void RHSFunction(Vector<T>& y, Vector<T>& yp)
		{
			double timedt = this->Info.time+this->Info.dt;
			this->Fyt = y;
			(*(this->Function))(timedt, y(0), this->Fyt(0));
			yp(0) = this->Y(0) - y(0) + (this->Info.dt)*(this->Fyt(0));
		};
		void LHSFunction(Vector<T>& z, Vector<T>& yp)
		{
			double timedt = this->Info.time+this->Info.dt;
			double dt = this->Info.dt;
			this->Fyt = z;
			(*(this->Jacobian))(timedt,this->Yk(0),z(0), this->Fyt(0));
			yp(0) = z(0)-dt*(this->Fyt(0));
			
		};
};


template<typename T>
class TrapezoidRule: public ODE_Subclass<T>
{
	public:
		TrapezoidRule(double t0,double t1,double dt):ODE_Subclass<T>(t0,t1,dt)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Fyt.SetDim(2);
			this->ExplicitMethod = 0;
		};
		TrapezoidRule(double t0,double t1,int nsteps):ODE_Subclass<T>(t0,t1,nsteps)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Fyt.SetDim(2);
			this->ExplicitMethod = 0;
		};
		void RHSFunction(Vector<T>& y, Vector<T>& yp)
		{
			double timedt = this->Info.time+this->Info.dt;
			this->Fyt(0) = y(0);
			this->Fyt(1) = y(0);
			(*(this->Function))(timedt, y(0), this->Fyt(0));
			(*(this->Function))(this->Info.time, this->Y(0),this->Fyt(1));
			yp(0) = this->Y(0) - y(0) + 0.5*(this->Info.dt)*(this->Fyt(0)+this->Fyt(1));
		};
		void LHSFunction(Vector<T>& z, Vector<T>& yp)
		{
			double timedt = this->Info.time+this->Info.dt;
			double dt = this->Info.dt;
			this->Fyt(0) = z(0);
			
			(*(this->Jacobian))(timedt,this->Yk(0),z(0), this->Fyt(0));
			yp(0) = z(0)-0.5*dt*(this->Fyt(0));
		};
};

template<typename T>
class RungeKutta4: public ODE_Subclass<T>
{
	public:
		RungeKutta4(double t0,double t1,double dt):ODE_Subclass<T>(t0,t1,dt)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Fyt.SetDim(4);
			this->ExplicitMethod = 1;
		};
		RungeKutta4(double t0,double t1,int nsteps):ODE_Subclass<T>(t0,t1,nsteps)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Fyt.SetDim(4);
			this->ExplicitMethod = 1;
		};
		void RHSFunction(Vector<T>& y, Vector<T>& yp)
		{
			this->Fyt(0) = y(0);
			this->Fyt(1) = y(0);
			this->Fyt(2) = y(0);
			this->Fyt(3) = y(0);
			double time = this->Info.time;
			double dt = this->Info.dt;
			double time_hdt = time+0.5*dt;
			double time_dt = time+dt;
			
			//stage 1
			(*(this->Function))(time, this->Y(0), this->Fyt(0));
			yp(0) = y(0)+0.5*dt*this->Fyt(0);
			//stage 2
			(*(this->Function))(time_hdt, yp(0), this->Fyt(1));
			yp(0) = y(0)+0.5*dt*this->Fyt(1);
			//stage 3
			(*(this->Function))(time_hdt, yp(0), this->Fyt(2));
			yp(0) = y(0)+dt*this->Fyt(2);
			//stage 4
			(*(this->Function))(time_dt, yp(0), this->Fyt(3));
			yp(0) = (1.0)/(6.0)*dt*(this->Fyt(0)+2.0*this->Fyt(1)+2.0*this->Fyt(2)+this->Fyt(3))+y(0);
			
		};
	private:
		Vector<T> Fyt;
};

#endif
