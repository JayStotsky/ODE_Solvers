#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include "Vector.H"
#include "ODE_Solver.H"


template<class T> FixedStep_ODE_Solver<T>::FixedStep_ODE_Solver(double tstart, double tend, int nsteps)
{
	this->T_start = tstart;
	this->T_start_set = 1;	
	
	this->T_end = tend;
	this->T_end_set = 1;
	
	Nsteps = nsteps;
	deltaT = (this->T_end-this->T_start)/(double)Nsteps;
	ready = 1;
};

template<class T> FixedStep_ODE_Solver<T>::FixedStep_ODE_Solver(double tstart, double tend, double dt)
{
	this->T_start = tstart;
	this->T_start_set = 1;	
	
	this->T_end = tend;
	this->T_end_set = 1;
	
	deltaT = dt;	
	Nsteps = floor((this->T_end-this->T_start)/dt);
	
	ready = 1;
};

template<class T> void FixedStep_ODE_Solver<T>::Initialize(double tstart, double tend, double dt)
{
	this->T_start = tstart;
	this->T_start_set = 1;
	
	this->T_end = tend;
	this->T_end_set = 1;
	
	deltaT = dt;
	Nsteps = floor((tend-tstart)/dt);
	
	ready = 1;
};

template<class T> void FixedStep_ODE_Solver<T>::Initialize(double tstart, double tend, int nsteps)
{
	this->T_start = tstart;
	this->T_start_set = 1;
	
	this->T_end = tend;
	this->T_end_set = 1;
	
	Nsteps = nsteps;
	deltaT = (tend-tstart)/(double)nsteps;
	
	ready =1;
};


template<class T> void FixedStep_ODE_Solver<T>::Solve()
{
	if (this->CheckStatus() && ready)
		for (int t = 0; t<Nsteps; t++)
			StepForward();
	this->SolverDone = 1;
}

template<class T> void FixedStep_ODE_Solver<T>::Solve(Vector<double>& Time, Vector<T>& Solution)
{
	if (this->CheckStatus() && ready)
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
			StepForward();
			
			Time(t+1) = this->time;
			Solution(t+1) = this->Y;
		}
	}
	this->SolverDone = 1;
}
