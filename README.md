Using ODE and PDE algorithms, we get numerical solutions to differential equations in thermodynamics.

# Bratu Equation
The bratu equation is a steady state model in thermodynamics. Its ODE is given by:  
<img src="https://render.githubusercontent.com/render/math?math=-\frac{d^2u}{dx}=\sigma*e^{u},u(0)=u(1)=0">
  
This was computationally solved by discretizing into a system of nonlinear equations using finite difference with <em>h</em> spacing.  
Then, using Newton's method this was solved.  
With Richardson Extrapolation, the solution is computed with higher accuracy for <em>h=0</em>.  
  
Stability and limits of this solution are analyzed in the code.  
  
# Diffusion
This is a simulation for the diffusion of heat distribution over a line. Temperature is then a function of linear position <em>x</em> and time <em>t</em>. Diffusion of a line is given by the PDE:  
<img src="https://render.githubusercontent.com/render/math?math=u_{t}=ku_{xx}">
Where <em>k</em> is a real number. The initial and boundary conditions I used were given by:
<img src="https://render.githubusercontent.com/render/math?math=u(x,0)=1,u(0,t)=u(1,t)=0">
  
This was solved by discretizing into a system of ODEs via finite difference with <em>h</em> spacing.  
Then, using Euler Forward, this ODE system was solved.
