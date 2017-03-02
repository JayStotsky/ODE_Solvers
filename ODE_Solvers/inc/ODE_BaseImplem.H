#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>

#include "Vector.H"
#include "Matrix.H"

template<typename T> ODE_Subclass<T>::ODE_Subclass(double t0,double t1,double dt)
{
	Info.dt = dt;
	Info.start_time = t0;
	Info.end_time = t1;
	Info.time = t0;
	Info.nsteps = floor(abs((t1-t0)/dt));
	Info.defined = 1;
	
};
template<typename T> ODE_Subclass<T>::ODE_Subclass(double t0,double t1,int nsteps)
{
	Info.nsteps = nsteps;
	Info.start_time = t0;
	Info.end_time = t1;
	Info.time = t0;
	Info.dt = (t1-t0)/(double)nsteps;
	Info.defined = 1;
};

template<typename T> void ODE_Subclass<T>::SetFunction(std::shared_ptr<std::function<void(double&,T&,T&)> > F)
{
	Function = F;
};
template<typename T> void ODE_Subclass<T>::SetJacobian(std::shared_ptr<std::function<void(double&,T&,T&,T&)> > J)
{
	Jacobian = J;
};

template<typename T> void ODE_Subclass<T>::SetInitialCondition(T& Y0)
{
	for (int i = 0; i < Y.GetDim(); i++)
		Y(i) = Y0;
	Yp = Y;
}


template<typename T> ODE_Superclass<T>::ODE_Superclass(std::shared_ptr<ODE_Subclass<T> > ts)
{
	TS = ts;
	NLS.SetODE(TS);
	ExplicitMethod = ts->ExplicitMethod;
};

	
template<typename T> void ODE_Superclass<T>::Solve(Vector<int>& t, Vector<Vector<T> >& y)
{
	if (TS->Info.defined)
	{
		store_time = &t;
		store_Y = &y;
		t(0) = TS->Info.start_time;
		y(0) = TS->Y;
		if(ExplicitMethod)
		{
			for (int k = 0; k < TS->Info.nsteps; k++)
			{
				SolveTS_Explicit();
				TS->Info.time = TS->Info.time+TS->Info.dt;
				
				Print(k);
			}
		}
		if(!ExplicitMethod)
		{
			for (int k = 0; k < TS->Info.nsteps; k++)
			{
				SolveTS_Implicit();
				TS->Info.time = TS->Info.time+TS->Info.dt;
				
				Print(k);
			}
		}	
	}
};


template<typename T> void ODE_Superclass<T>::SolveTS_Implicit()
{
	TS->Yp = TS->Y;
	NLS.Solve(TS->Yp);
	TS->Y = TS->Yp;
};


template<typename T> void ODE_Superclass<T>::SolveTS_Explicit()
{
	TS->RHSFunction(TS->Y,TS->Yp);
	TS->Y = TS->Yp;
};

template<typename T> void ODE_Superclass<T>::Print(int k)
{
	for (int i = 0; i < store_time->GetDim(); i++)
		if (k+1==(*store_time)(i))
			(*store_Y)(i)=TS->Y;
}
