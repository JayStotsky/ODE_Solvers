#ifndef _FIXEDSTEPMETHODS_H_
#define _FIXEDSTEPMETHODS_H_


#include "FixedStep_ODE_Solvers.H"
#include "FunctionClass.H"
#include "Vector.H"

template<class T> class FixedStepForwardEuler: public FixedStep_ODE_Solver<T>{
	public:		
		FixedStepForwardEuler(){};
		FixedStepForwardEuler(double tstart, double tend, int nsteps): FixedStep_ODE_Solver<T>(tstart, tend, nsteps){}; //call constructor in base class
		FixedStepForwardEuler(double tstart, double tend, double dt): FixedStep_ODE_Solver<T>(tstart, tend, dt){}; //call constructor in base class
		
		void StepForward()
		{
			(*(this->function))(this->time,this->Y,this->Fty);
			this->Y = this->Y+(this->deltaT*(this->Fty));
			this->time = this->time+this->deltaT;
		};
	private:
};


template<class T> class FixedStepRK4: public FixedStep_ODE_Solver<T>{
	public:		
		FixedStepRK4(){};
		FixedStepRK4(double tstart, double tend, int nsteps): FixedStep_ODE_Solver<T>(tstart, tend, nsteps){}; //call constructor in base class
		FixedStepRK4(double tstart, double tend, double dt): FixedStep_ODE_Solver<T>(tstart, tend, dt){}; //call constructor in base class
		
		void StepForward()
		{
			k1=this->Y;
			k2=this->Y;
			k3=this->Y;
			k4=this->Y;
			double deltaT = this->deltaT;
			
			(*(this->function))(this->time,this->Y,k1);
					
			f1 = this->Y+0.5*deltaT*k1;
			t1 = this->time+0.5*deltaT;
			(*(this->function))(t1,f1,k2);
			
			t2 = this->time+0.5*deltaT;
			f2 = this->Y+0.5*deltaT*k2;
			(*(this->function))(t2,f2,k3);
			
			t3 = this->time+deltaT;
			f3 = this->Y+deltaT*k3;
			(*(this->function))(t3,f3,k4);
			
			this->Fty = (1.0/6.0)*(k1+2*k2+2*k3+k4);
			this->Y = this->Y+(deltaT*this->Fty);
			this->time = this->time+deltaT;
		};
	private:
		T k1, k2, k3, k4;
		T f1, f2, f3;
		double t1, t2, t3;
		
};
	
template<class T> class ImplicitFixedStep: public FixedStep_ODE_Solver<T>{
	public:
		ImplicitFixedStep(double tstart, double tend, int nsteps) : FixedStep_ODE_Solver<T>(tstart, tend, nsteps){};
		ImplicitFixedStep(double tstart, double tend, double dt) : FixedStep_ODE_Solver<T>(tstart, tend, dt){};
		void SetJacobianFunction(std::shared_ptr<std::function<void(double&,T&,T&,T&)> > J)
		{
			Jacobian = J;
			JacobianSet = 1;
		}
		
		void ConstructJacobian(T& yk){Yk = yk;};	
		virtual void EvalJacobian(T& in,T& out){};
		virtual void EvalRightHandSide(T& in,T& out){};
		void StepForward();
		void SetNewtonSolver(std::shared_ptr<ODENewton<T> > NS)
		{ N = NS;};
		
	protected:
		std::shared_ptr<std::function<void(double&,T&,T&,T&)> > Jacobian;
		bool JacobianSet=0;
		ODENewton<T> N;
		T Yk;

};
	
	
template<class T> class BackwardEuler: public ImplicitFixedStep<T>{
	public:
		BackwardEuler(double tstart, double tend, int nsteps) : ImplicitFixedStep<T>(tstart, tend, nsteps){};
		BackwardEuler(double tstart, double tend, double dt) : ImplicitFixedStep<T>(tstart, tend, dt){};
		void EvalJacobian(T& in,T& out)
		{
			double time = this->time+this->deltaT;
			(*(this->Jacobian))(time, this->Yk, in, out);
			out = in-this->deltaT*out;
		};
		void EvalRightHandSide(T& in,T& out)
		{
			double time = this->time+this->deltaT;
			(*(this->function))(time, in, out);
			out = this->Y+this->deltaT*out-in;
		};
	private:
};


template<class T> class Trapezoid: public ImplicitFixedStep<T>{
	public:
		Trapezoid(double tstart, double tend, int nsteps) : ImplicitFixedStep<T>(tstart, tend, nsteps){};
		Trapezoid(double tstart, double tend, double dt) : ImplicitFixedStep<T>(tstart, tend, dt){};

		void EvalJacobian(T& in,T& out)
		{
			double time = this->time+this->deltaT;
			(*(this->Jacobian))(time, this->Yk, in, out);
			out = in-0.5*this->deltaT*out;
		};
		void EvalRightHandSide(T& in,T& out)
		{
			T fY;
			double time = this->time+this->deltaT;
			(*(this->function))(this->time, in, out);
			(*(this->function))(time, this->Y, fY);
			out = this->Y+0.5*(this->deltaT*out+this->deltaT*fY)-in;
		};
	private:
};

	
#include "ImplicitMultistepImplem.H"

#endif 
