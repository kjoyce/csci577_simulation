from scipy.integrate import odeint  # for integrate.odeint
from pylab import *  # for plotting commands
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib.mlab import find
debug_here = Tracer()
 
# Notice that parameters are global.
GM = 4*pi**2
Gm1 = .001*GM
Gm2 = .04*GM
years = 100.
num_samples = 1000
 
def two_body(state,t):  #OUCH! The signature is reversed for odeint!
  x1 = array([state[0],state[1]])  # body 1 position
  v1 = array([state[2],state[3]])  # body 1 velocity
  x2 = array([state[4],state[5]])
  v2 = array([state[6],state[7]])

  r1  = norm(x1)
  r2  = norm(x2)
  r21 = norm(x2 - x1)

  a1 = -GM/r1**3*x1 + Gm2/r21**3*(x2-x1)
  a2 = -GM/r2**3*x2 + Gm1/r21**3*(x1-x2)
  derivatives = array([v1,a1,v2,a2]).T.flatten(1)
  return derivatives
 
times = linspace(0.0,years,num_samples)
xinit = array([2.52,0,0,sqrt(GM/2.52),5.24,0,0,sqrt(GM/5.24)])  # initial values
x = odeint(two_body,xinit,times)
x = x.T

r1 = sqrt(x[0]**2 + x[1]**2)
r2 = sqrt(x[4]**2 + x[5]**2)
r21 = sqrt( (x[4]-x[0])**2 + (x[5]-x[1])**2)
E = .5*(Gm1*x[2]**2 + Gm2*x[3]**2) - GM*( Gm1/r1 + Gm2/r2 + Gm1*Gm2/r21 ) 
L = Gm1*(x[0]*x[3] - x[1]*x[2]) + Gm2*(x[4]*x[7] - x[5]*x[6])

deltaE = (E[0]-E)/E[0]
deltaL = (L[0]-L)/L[0]

subplot(311)
plot(x[0],x[1])
plot(x[4],x[5])
xlim((-20,20))
grid()
title('Planetary Orbit')
xlabel('Horizontal Distance (AU)')
ylabel('Vertical Distance (AU)')

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

show()
