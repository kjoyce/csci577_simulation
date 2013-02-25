from matplotlib.pyplot import *
import Nbody 
from scipy.integrate import odeint
from numpy import *

G = .5
reader = Nbody.InitialConditions()
sys_1 = reader.read_file('euler_init.txt')

force = Nbody.Force(sys_1,G)

t = linspace(0.,100.,1000)
x = odeint(force,sys_1.init,t)
print shape(x)
plot(x[:,0],x[:,1])
plot(x[:,2],x[:,3])
plot(x[:,4],x[:,5])
show()
