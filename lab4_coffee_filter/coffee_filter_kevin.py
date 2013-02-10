from scipy import *
from numpy import *
from numpy.linalg import inv
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib.mlab import find
from ODE_integrator import RungeKutta 
debug_here = Tracer()
data = array([[.2055,.2302,.2550,.2797,.3045,.3292,.3539,.3786,.4033,.4280,
              .4526,.4773,.5020,.5266,.5513,.5759,.6005,.6252,.6498,.6744,
              .6990,.7236,.7482,.7728,.7974,.8220,.8466],
              [.4188,.4164,.4128,.4082,.4026,.3958,.3878,.3802,.3708,.3609,
               .3505,.3400,.3297,.3181,.3051,.2913,.2788,.2667,.2497,.2337,
               .2175,.2008,.1846,.1696,.1566,.1393,.1263]])

t = data[0]
y = data[1]

def numeric_derivative(t,y):
  tt = (t[0:-1] + t[1:])/2
  yy = (y[1:] - y[0:-1])/(t[1:] - t[0:-1])
  return (tt,yy)

(tt,yy) = numeric_derivative(t,y)
(ttt,yyy) = numeric_derivative(tt,yy)
subplot(221)
plot(t,y,'ro-')
subplot(222)
plot(tt,yy,'go-')
subplot(223)
plot(ttt,yyy,'ko-')
subplot(224)
plot(yy[1:],yyy,'ko')

# solve linear least squares problem
vt = (1/dot(yy,yy)*dot(yy[1:],yyy/9.8 + 1))**(-1) 
#vt = -.45

def f1(t,y):
  return array([y[1], -9.8*(1-y[1]/vt)])

def f2(t,y):
  return array([y[1], -9.8*(1-(y[1]/-vt**2)**2)])

# This is a new model that has both v and v**2,  Not much better
# solve least squares
A = array([yy[1:],yy[1:]**2])
[k1, k2] = dot(inv(dot(A,A.T)),dot(A,yyy+9.8))

def f3(t,y):
  return array([y[1], -9.8 + k1*y[1] + k2*y[1]**2])


# This is a nonlinear model that has the parameter a
# F = F_g + k v ^ a
# Something is fishy with my least squares parameter selection
# though
A=vstack([ones(size(yyy)),log(-yy[1:])])
[k,a] = dot(inv(dot(A,A.T)),dot(A,log(yyy+9.8)))
k = exp(k)
def f4(t,y):
  return array([y[1], -9.8 + -k*(-y[1])**a])

for f in (f1,f2,f3,f4):
  integrator = RungeKutta(f,[y[0],0],data[0][0],data[0][-1],.001)
  (mt,my) = integrator.integrate()
  idx = find(my[0] > 0)
  mt = mt[idx]
  my = my[:,idx]
  figure()
  plot(t,y,'.')
  plot(mt,my[0])
  show()


# use fmin to find parameters
