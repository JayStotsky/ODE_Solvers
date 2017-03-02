#ifndef _DISCRETE_GALERKIN_TS_H_
#define _DISCRETE_GALERKIN_TS_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"
#include "ODE_Quadrature.H"
#include "ODE_Base.H"


template<typename T>
class DiscreteGalerkin0:public ODE_Subclass<T>
{
	public:
		DiscreteGalerkin0(double t0,double t1,double dt):ODE_Subclass<T>(t0,t1,dt)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Yk.SetDim(1);
			this->Fyt.SetDim(1);
			this->ExplicitMethod = 0;
		};
		DiscreteGalerkin0(double t0,double t1,int nsteps):ODE_Subclass<T>(t0,t1,nsteps)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Yk.SetDim(1);
			this->Fyt.SetDim(1);
			this->ExplicitMethod = 0;			
		};
		
		void RHSFunction(Vector<T>& y, Vector<T>& yp)
		{
			I = y(0)*0.0;
			this->Yk(0) = y(0);
			double tm = this->Info.time;
			double tp = this->Info.time + this->Info.dt;
			Q.Integrate(tm, tp, &DiscreteGalerkin0::Ft_Yt, I, this); 
			yp(0) = this->Y(0)+I-y(0) ;
			
		};
		void LHSFunction(Vector<T>& y, Vector<T>& yp)
		{
			I = y(0)*0.0;
			double tm = this->Info.time;
			double tp = this->Info.time + this->Info.dt;
			this->Yp(0) = y(0);
			Q.Integrate(tm, tp, &DiscreteGalerkin0::Jt_Yt, I, this);
			yp(0) = y(0) - I;
			
		};
		
		void Ft_Yt(double& t, T& F)
		{
			(*(this->Function))(t,this->Yk(0),F);
		};
		
		void Jt_Yt(double& t, T& F)
		{
			(*(this->Jacobian))(t, this->Yk(0), this->Yp(0), F);
		};
		
		void SetQuadrature(std::shared_ptr<ODE_Quadrature_Base<T,DiscreteGalerkin0<T> > > q)
		{Q = *q;};
	
	
	private:
		ODE_Quadrature_Base<T,DiscreteGalerkin0<T> > Q;
		T I;
	
};

class DiscreteGalerkin1:public ODE_Subclass<T>
{
	public:
		DiscreteGalerkin1(double t0,double t1,double dt):ODE_Subclass<T>(t0,t1,dt)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Yk.SetDim(1);
			this->Fyt.SetDim(1);
			this->ExplicitMethod = 0;
		};
		DiscreteGalerkin1(double t0,double t1,int nsteps):ODE_Subclass<T>(t0,t1,nsteps)
		{
			this->Y.SetDim(1);
			this->Yp.SetDim(1);
			this->Yk.SetDim(1);
			this->Fyt.SetDim(1);
			this->ExplicitMethod = 0;			
		};
		
		void RHSFunction(Vector<T>& y, Vector<T>& yp)
		{
			I = y(0)*0.0;
			this->Yk(0) = y(0);
			double tm = this->Info.time;
			double tp = this->Info.time + this->Info.dt;
			Q.Integrate(tm, tp, &DiscreteGalerkin1::Ft_Yt, I, this); 
			yp(0) = this->Y(0)+I-y(0) ;
			
		};
		void LHSFunction(Vector<T>& y, Vector<T>& yp)
		{
			I = y(0)*0.0;
			double tm = this->Info.time;
			double tp = this->Info.time + this->Info.dt;
			this->Yp(0) = y(0);
			Q.Integrate(tm, tp, &DiscreteGalerkin1::Jt_Yt, I, this);
			yp(0) = y(0) - I;
			
		};
		
		void Ft_Yt(double& t, T& F)
		{
			(*(this->Function))(t,this->Yk(0),F);
		};
		
		void Jt_Yt(double& t, T& F)
		{
			(*(this->Jacobian))(t, this->Yk(0), this->Yp(0), F);
		};
		
		void SetQuadrature(std::shared_ptr<ODE_Quadrature_Base<T,DiscreteGalerkin1<T> > > q)
		{Q = *q;};
	
	
	private:
		ODE_Quadrature_Base<T,DiscreteGalerkin1<T> > Q;
		T I;
	
};


#endif
