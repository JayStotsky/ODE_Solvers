#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include "Vector.H"
#include "ODE_Solver.H"


template<class T, class S> RK4<T,S>::RK4(double tstart, double tend, int nsteps)
{
	this->T_start = tstart;
	this->T_start_set = 1;	
	
	if (tend>tstart)
	{
		this->T_end = tend;
		this->T_end_set = 1;
	}
	
	if (nsteps>0)
		Nsteps = nsteps;
		
	if (this->T_end_set && Nsteps>0)
		deltaT = (this->T_end-this->T_start)/(double)Nsteps;
		
	this->Y.SetDim(Nsteps+1);
	this->time.SetDim(Nsteps+1);
};

template<class T, class S> RK4<T,S>::RK4(double tstart, double tend, double dt)
{
	this->T_start = tstart;
	this->T_start_set = 1;	
	
	if (tend>tstart)
	{
		this->T_end = tend;
		this->T_end_set = 1;
	}
	
	if (dt>0)
		deltaT = dt;
		
	if (this->T_end_set && deltaT>0)
		Nsteps = floor((this->T_end-this->T_start)/dt);
		
	this->Y.SetDim(Nsteps+1);
	this->time.SetDim(Nsteps+1);
};

template<class T, class S> void RK4<T,S>::Initialize(double tstart, double tend, int nsteps)
{
	this->T_start = tstart;
	this->T_start_set = 1;
	
	if (tend>tstart)
	{
		this->T_end = tend;
		this->T_end_set = 1;
	}
	
	if (nsteps>0)
		Nsteps = nsteps;
		
	this->Y.SetDim(Nsteps+1);
	this->time.SetDim(Nsteps+1);
};

template<class T, class S> void RK4<T,S>::Solve()
{
	if (this->CheckStatus())
	{
		T k1, k2, k3, k4;
		T f1, f2, f3;
		T t1, t2, t3;
		this->time(0) = this->T_start;
		for (int t = 0; t<Nsteps; t++)
		{
			this->time(t+1) = this->time(t)+deltaT;
			
			this->Function->eval(this->time(t),this->Y(t),k1);
			
			f1 = this->Y(t)+0.5*deltaT*k1;
			t1 = this->time(t)+0.5*deltaT;
			this->Function->eval(t1,f1,k2);
			
			t2 = this->time(t)+0.5*deltaT;
			f2 = this->Y(t)+0.5*deltaT*k2;
			this->Function->eval(t2,f2,k3);
			
			t3 = this->time(t)+deltaT;
			f3 = this->Y(t)+0.5*deltaT*k3;
			this->Function->eval(t3,f3,k4);
			
			this->Fty = (1.0/6.0)*(k1+2*k2+2*k3+k4);
			this->Y(t+1) = this->Y(t)+(deltaT*this->Fty);
		}
	}
	this->SolverDone = 1;
}
