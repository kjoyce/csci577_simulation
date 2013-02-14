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
 
times = linspace(0.0,20.,300)
yinit = array([2.52,0,0,sqrt(GM/2.52),5.24,0,0,sqrt(GM/5.24)])  # initial values
y = odeint(two_body,yinit,times)
y = y.T
#E = .5*Gm1*sum(y[2:3]**2) - Gm1/r1 - Gm2/r2 - 

subplot(311)
plot(y[0],y[1])
plot(y[4],y[5])
xlim((-20,20))
grid()
title('Planetary Orbit')
xlabel('Horizontal distance (AU)')
ylabel('Vertical Distance (AU)')

subplot(312)
grid()

subplot(313)
grid()

show()
