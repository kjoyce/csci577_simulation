from matplotlib.pyplot import *
import Nbody 
from scipy.integrate import odeint
from numpy import *

G = .5
files = ('euler_init.txt','montgomery_init.txt','lagrange_init.txt')
for fname in files:
  reader = Nbody.InitialConditions()
  sys_1 = reader.read_file(fname)
  force = Nbody.Force(sys_1,G)
  t = linspace(0.,10.,100)
  x = odeint(force,sys_1.init,t)
  plot(x[:,0],x[:,1])
  plot(x[:,2],x[:,3])
  plot(x[:,4],x[:,5])
  show()
