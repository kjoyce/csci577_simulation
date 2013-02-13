from scipy.integrate import odeint  # for integrate.odeint
from pylab import *  # for plotting commands
 
# Notice that parameters are global.
GM = 4*pi**2
m1 = .001*GM
m2 = .04*GM
 
def two_body(state,t):  #OUCH! The signature is reversed for odeint!
  (x1,y1) = (state[0],state[1])  # body 1 position
  (v1,w1) = (state[2],state[3])  # body 1 velocity
  (x2,y2) = (state[4],state[5])
  (v2,w2) = (state[6],state[7])

  ax1 = GM/(x1**2 + y1**2)**(1.5)*(x2 - x1)
  ay1 = GM/(x1**2 + y1**2)**(1.5)*(y2 - y1)
  return 
 
#times = linspace(0.0,5.*2.*pi*sqrt(k/m),1000)
#yinit = array([1.,0.])  # initial values
#y = odeint(sho,yinit,times)
# 
## Plot the results
#clf()
#plot(times,y[:,0]) # y[:,0] is the first column of y, the postitions
#xlabel('Times (s)')
#ylabel('Positions (m)')
#title('Simple Harmonic Motion')
#grid()
#show()
