import Nbody 
from scipy.integrate import odeint
from numpy import *

G = 6.67*10**(-11)
reader = Nbody.InitialConditions()
sys_1 = reader.read_file('euler_init.txt')

force = Nbody.Force(sys_1,G)

t = linspace(0.,100.,1000)
x = odeint(force,sys_1.init,t)
