#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include "Vector.H"
#include "ODE_Solver.H"


template<class T, class S> ForwardEuler<T,S>::ForwardEuler(double tstart, double tend, int nsteps)
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

template<class T, class S> ForwardEuler<T,S>::ForwardEuler(double tstart, double tend, double dt)
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

template<class T, class S> void ForwardEuler<T,S>::Initialize(double tstart, double tend, int nsteps)
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

template<class T, class S> void ForwardEuler<T,S>::Solve()
{
	if (this->CheckStatus())
	{	
		this->time(0) = this->T_start;
		for (int t = 0; t<Nsteps; t++)
		{
			this->time(t+1) = this->time(t)+deltaT;
			this->Function->eval(this->time(t),this->Y(t),this->Fty);
			this->Y(t+1) = this->Y(t)+(deltaT*this->Fty);
			if (typeid(this->Y(t+1))==typeid(Vector<double>))
			{
				Vector<double> V = this->Y(t+1);
				printf("%f %f %d\n", V(0),V(1),V.GetDim());
			}
		}
	}
	this->SolverDone = 1;
}
