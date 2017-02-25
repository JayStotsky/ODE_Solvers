Simple solvers for ODE's. Goal is to have solvers that rely on template classes to maintain generality.
The ODE 
y'=f(t,y)
y(0)=y0  
is solved in each case. The function, f(t,y) should be made from a class with a member function called "eval". The class-type is a template argument of the ODE solver. 

1) Explicit Multistep Methods (constant timestep)
2) Explicit Runge-Kutta Methods (constant timestep)

To be added:
Discrete Galerkin time step routines
Implicit multistep
Implicit Runge-Kutta
adaptivity
