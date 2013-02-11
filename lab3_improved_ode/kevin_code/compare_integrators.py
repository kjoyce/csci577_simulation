from scipy import *
from numpy import *
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib.mlab import find
debug_here = Tracer()
from ODE_integrator import RungeKutta, Euler, EulerRichardson, PredictorCorrector

def oscillator_system(y,k,m):
  return array([y[1], -k*y[0]])

def oscillator_analytic(t,k,x0,m):
  return array([x0*cos(sqrt(k/m)*t), -x0*sqrt(k/m)*sin(sqrt(k/m)*t)]).T

def oscillator_energy(x,v,k,m):
  return .5*k*x**2 + .5*m*v**2

def oscillator_error(x,v,k,m):
  e = oscillator_energy(x,v,k,m)
  return abs((e[0] - e[-1])/e[0])

# The purpose of this class is to make indexing a little more readable
# I.e. rather than referring to results[0][1][4] or something like that, you can
# >tester = TestIntegrator(RungeKutta,.01)
# >tester.iterate_test()
# >print tester.n # this says how many iterations it took
# >print tester.dts[-1] # this is the last dt needed
class TestIntegrator:
  def __init__(self,method,tolerance=.01):
    self.method = method  
    self.tolerance = tolerance
    self.k = 1		  # spring constant	
    self.m = 1		  # mass
    (self.a,self.b) = (0.,10*pi) # five oscillations
    self.y0 = array([1,0]) # initial position
    self.f = lambda t,y:oscillator_system(y,self.k,self.m) # make system jive with integrator
    self.n = 0 # initialize count for integrating 
    self.dts = [] #initialize dts
    self.integrators = [] # initialize list of ODE_integrators
    self.errors = [] # initialize energy errors

  def iterate_test(self,next_dt):
    while True:
      self.n += 1
      dt = next_dt(self.n) #next_dt needs to be a decreasing function of n
      integrator = self.method(self.f,self.y0,self.a,self.b,dt)
      (t,y) = integrator.integrate()
      self.integrators.append(integrator)
      self.dts.append(dt)
      err = oscillator_error(y[0],y[1],self.k,self.m)
      self.errors.append(err)
      print "dt = {}, error = {}".format(dt,err)
#      if self.n>2 and self.errors[-1] > self.errors[-2]:
#	print "Non-Decreasing Errors!"
#	break
      if self.errors[-1] < self.tolerance:
	break
    return

  # you can pass some function of dt
  def plot_errors(self,f=lambda dt:dt,upper=Inf):
    dts = array(self.dts)
    errs = array(self.errors)
    idx = find(dts<upper)
    p = plot(f(array(dts[idx])),errs[idx],'o')
    return p

### the following function change the density of dts
# this decreases by .5, .01, .005,...
def exp_decrease(n):
  n = n+2
  return (4.*((n-1)%2)+1)/10**(n/2) 

# this is 10, 
def linear_decrease(n):
  return 1./(n*10.)*10

# this is the logarithm of linear_decrease
def log_decrease(n):
    return 1./log(1./linear_decrease(n+1))

# Usage example
tolerance = .01
print "----- Runge-Kutta --------"
rungekutta = TestIntegrator(RungeKutta,tolerance)
rungekutta.iterate_test(log_decrease)

print "----- Euler-Richardson --------"
eulerrich = TestIntegrator(EulerRichardson,tolerance)
eulerrich.iterate_test(linear_decrease)

print "----- Predictor Corrector --------"
pc = TestIntegrator(PredictorCorrector,tolerance)
pc.iterate_test(linear_decrease)

# Be careful with this... it takes a while at 10^-6
print "----- Naive Euler --------"
euler = TestIntegrator(Euler,tolerance)
euler.iterate_test(exp_decrease)

