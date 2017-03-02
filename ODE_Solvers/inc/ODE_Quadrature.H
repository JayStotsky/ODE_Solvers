#ifndef _ODE_QUADRATURE_H_
#define _ODE_QUADRATURE_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"
#include "ODE_Base.H"

template<typename T> class DiscreteGalerkin0;

template<typename T, class S>
class ODE_Quadrature_Base
{
	public:
		ODE_Quadrature_Base(){};
		ODE_Quadrature_Base(double L,double U)
		{
			lower_bound = L;
			upper_bound = U;
		};
		void SetBounds(double L,double U)
		{
			lower_bound = L;
			upper_bound = U;
		};
		void SetPointsWeights(Vector<double>& P, Vector<double>& W)
		{
			points = P;
			weights = W;
		};
		
		void Integrate(void(S::*F)(double&,T&), T& I, S* DG)
		{
			T f = I;
			double t;
			double range = upper_bound - lower_bound;
			for (int i = 0; i < points.GetDim(); i++)
			{
				t = range*points(i) + lower_bound;
				(DG->*F)(t,f);
				I += f * weights(i) * range;
			}
			
		};
		
		void Integrate(double L, double U,void(S::*F)(double&,T&), T& I, S* DG)
		{
			T f = I; //initialize f
			double t;
			double range = U - L;
			for (int i = 0; i < points.GetDim(); i++)
			{
				t = range*points(i) + L;
				(DG->*F)(t,f);
				I += f * weights(i) * range;
			}
			
		};
		
	protected:
		Vector<double> points;
		Vector<double> weights;
		double lower_bound;
		double upper_bound;
};

template<typename T, class S>
class MidpointRuleQuad:public ODE_Quadrature_Base<T,S>
{
	public:
		MidpointRuleQuad()
		{
			Vector<double> points(1);
			Vector<double> weights(1);
			points(0) = 0.5;
			weights(0) = 1.0;
			this->SetPointsWeights(points,weights);
		};
		MidpointRuleQuad(double L,double U):ODE_Quadrature_Base<T,S>(L,U)
		{
			Vector<double> points(1);
			Vector<double> weights(1);
			points(0) = 0.5;
			weights(0) = 1.0;
			this->SetPointsWeights(points,weights);
		};
};


template<typename T, class S>
class TrapezoidRuleQuad:public ODE_Quadrature_Base<T,S>
{
	public:
		TrapezoidRuleQuad()
		{
			Vector<double> points(2);
			Vector<double> weights(2);
			points(0) = 0.0;
			points(1) = 1.0;
			weights(0) = 0.5;
			weights(1) = 0.5;
			this->SetPointsWeights(points,weights);
		};
		TrapezoidRuleQuad(double L,double U): ODE_Quadrature_Base<T,S>(L,U)
		{
			Vector<double> points(2);
			Vector<double> weights(2);
			points(0) = 0.0;
			points(1) = 1.0;
			weights(0) = 0.5;
			weights(1) = 0.5;
			this->SetPointsWeights(points,weights);
		};
};

template<typename T,class S>
class TrapezoidCompositeQuad:public ODE_Quadrature_Base<T,S>
{
	public:
		TrapezoidCompositeQuad(){};
		TrapezoidCompositeQuad(double L,double U, int n): ODE_Quadrature_Base<T,S>(L,U)
		{
			Vector<double> points(n);
			Vector<double> weights(n);
			
			for (int i = 0; i < n; i++)
			{
				points(i) = (double)i/(double)(n-1);
				weights(i) = 1.0/(double)n;
			}
			this->SetPointsWeights(points,weights);
		};
};



#endif
