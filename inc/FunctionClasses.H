#ifndef _FUNCTIONCLASSES_H_
#define _FUNCTIONCLASSES_H_

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"


class VectorFunction{
	public:
		virtual void eval(double &, Vector<double>&, Vector<double>&){};
	private:
};

class LinearSystem: public VectorFunction{
	public:
		
		void eval(double& t, Vector<double>& y, Vector<double>& f)
		{
			f = M*y;
		};
		void set(Matrix<double>& A)
		{M = A;};
	private:
		Matrix<double> M;
};

class NonlinearSystem: public VectorFunction{
	public:
		void eval(double &t, Vector<double>& y, Vector<double>& f)
		{
			(*NLS)(t,y,f);
		};
		
		void set(std::shared_ptr<std::function<void(double&,Vector<double>&,Vector<double>&)> > nls)
		{
			NLS = nls;
		}
	private:
		std::shared_ptr<std::function<void(double&,Vector<double>&,Vector<double>&)> > NLS;
};



template<class T> class FunctionWrapper
{
	public:
		void operator()(double& t, T& in, T& out)
		{(*F)(t,in,out);};
	private:
		std::shared_ptr<std::function<void(double&, T&, T&)> > F;
		
};



template<class T> class JacobianFunction
{
	public:
		void eval (T& in, T& out)
		{
			(*J)(t,y,in,out);
			out = in-deltaT*out ;
		};
		
		void setY(T& Y)
		{
			y = Y;
		};
		
		void setT(double& tt, double& dt)
		{
			t = tt;
			deltaT = dt;
		};
		
		void SetFunction(std::shared_ptr<std::function<void(double&, T&, T&, T&)> > F )
		{ J = F;};
		
	protected:
		std::shared_ptr<std::function<void(double&, T&, T&, T&)> > J;
		T y;
		double t, deltaT;
};


template<class T> class Function
{
	public:
		virtual void eval (T&, T&){};
		void setT(double& tt)
		{
			t = tt;
		};
		
		void SetFunction(std::shared_ptr<std::function<void(double&, T&, T&)> > F ){ J = F;};
		
	protected:
		std::shared_ptr<std::function<void(double&, T&, T&)> > J;
		double t;
};

template<class T> class RHSFunction: public Function<T>
{
	public:
		void eval( T& in, T& out)
		{
			(*(this->J))(this->t,in,out);

			out = y+deltaT*out-in;
		};
		void set(T& Y, double& dt, double& tt)
		{
			y = Y;
			deltaT = dt;
			this->t= tt;
		};
	private:
		T y;
		double deltaT;
};	


#endif

