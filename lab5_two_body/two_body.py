from scipy.integrate import odeint  # for integrate.odeint
from pylab import *  # for plotting commands
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib.mlab import find
from ODE_integrator import RungeKutta
from two_body_params import *
debug_here = Tracer()
def two_body(state,t):  #OUCH! The signature is reversed for odeint!
  x1 = array([state[0],state[1]])  # body 1 position
  v1 = array([state[2],state[3]])  # body 1 velocity
  x2 = array([state[4],state[5]])
  v2 = array([state[6],state[7]])

  r1  = norm(x1)
  r2  = norm(x2)
  r21 = norm(x2 - x1)

  a1 = -GM/r1**3*x1 + m2*GM/r21**3*(x2-x1)
  a2 = -GM/r2**3*x2 + m1*GM/r21**3*(x1-x2)
  derivatives = array([v1,a1,v2,a2]).T.flatten(1)
  return derivatives
 
scipy_result = odeint(two_body,xinit,times)
scipy_result = scipy_result.T

atoll=10**-13
rtoll=10**-14
high_prec_result = odeint(two_body,xinit,times,atol=atoll,rtol=rtoll)
high_prec_result = high_prec_result.T

######## Do RK4 HERE ############
f = lambda t,x: two_body(x,t)
rk4 = RungeKutta(f,xinit,0.,years,years/num_samples)
(t,rk4_result) = rk4.integrate()
rk4_result = rk4_result[:,:-2]
#################################

iterates = ((scipy_result,"Scipy odeint default tolerances"),
	    (rk4_result,"Runge Kutta 4th Order"),
	    (high_prec_result,"odeint rtol={} atol={}".format(rtoll,atoll)))

for x,t in iterates:
  figure()
  r1 = sqrt(x[0]**2 + x[1]**2)
  r2 = sqrt(x[4]**2 + x[5]**2)
  r21 = sqrt( (x[4]-x[0])**2 + (x[5]-x[1])**2) 
  E = .5*(m1*(x[2]**2 + x[3]**2) + m2*(x[6]**2+x[7]**2)) - GM*( m1/r1 + m2/r2 + m1*m2/r21 ) 
  deltaE = (E-E[0])/E[0]

  L = m1*(x[0]*x[3] - x[1]*x[2]) + m2*(x[4]*x[7] - x[5]*x[6])
  deltaL = (L-L[0])/L[0]

  subplot(311)
  plot(x[0],x[1])
  plot(x[4],x[5])
  xlim((-20,20))
  grid()
  title(t)
  xlabel('Horizontal Distance (AU)')
  ylabel('Vertical Distance (AU)')
  gca().axis('equal')

  subplot(312)
  plot(times,deltaE)
  xlabel('Time in years')
  ylabel('% $\Delta E /M$')
  grid()

  subplot(313)
  plot(times,deltaL)
  xlabel('Time in years')
  ylabel('% $\Delta L / M$')
  grid()

  subplots_adjust(wspace=.4)  # Note this makes space
  subplots_adjust(hspace=.4)  # Note this makes space
  show()


