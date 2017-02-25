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
};

template<class T, class S> void ForwardEuler<T,S>::Solve()
{
	if (this->CheckStatus())
	{	
		for (int t = 0; t<Nsteps; t++)
		{
			this->Function->eval(this->time,this->Y,this->Fty);
			this->Y = this->Y+(deltaT*(this->Fty));
			this->time = this->time+deltaT;
		}
	}
	this->SolverDone = 1;
}

template<class T, class S> void ForwardEuler<T,S>::Solve(Vector<double>& Time, Vector<T>& Solution)
{
	if (this->CheckStatus())
	{	
		for (int t = 0; t<Nsteps; t++)
		{	
			if (Time.GetDim()==0 || Solution.GetDim()==0)
			{
				Time.SetDim(Nsteps+1);
				Solution.SetDim(Nsteps+1);
			}
			
			if (t==0)
			{
				Time(0) = this->T_start;
				Solution(0) = this->Y;
			}
			
			this->Function->eval(this->time,this->Y,this->Fty);
			this->Y = this->Y+(deltaT*(this->Fty));
			this->time = this->time+deltaT;
			
			
			Time(t+1) = this->time;
			Solution(t+1) = this->Y;
		}
	}
	this->SolverDone = 1;
}
