"""
Simulation of Ice-Sheet temperature distribution
""" 
from dolfin import *
import numpy as np
import pylab as py

spy = 31556926.
# Parameters with meaningful change 
zb = -300.
v_ice_rate = -.2
slope = .1E-5
dth_dx = 2E-4
qgeo = 30E-3
n_mesh = 20
total_time = 1000*spy
dt = spy

# Constants 
zs = 0.
g = 9.81
spy = 31556926.
rho = 911.
cp = 2009.
beta = 9.8E-8
thetapmp = beta*rho*g*(zs - zb)
k = 2.1

# Surface temperature
def surface_temp(t):
  return -10 + 5*np.sin(2*pi*t/spy)

# Expressions
#sigma = lambda z: (z-zb)/(zs-zb)
#dsigma_dz = (zs-zb)**(-1)
#u = lambda z: 20-100*sigma(z)**4
#du_dz = lambda z: -400.*sigma(z)*dsigma_dz
#w = lambda z: v_ice*sigma(z)**4
# phi = lambda z: -rho*g*(zs - z)*du_dz(z)*slope
sigma_exp = '(x[0] - zb)/(zs - zb)'
sigma = Expression(sigma_exp,zb=zb,zs=zs)
u = Expression('20-100*pow('+sigma_exp+',4)',zb=zb,zs=zs)
w = Expression(('v_ice_rate*pow('+sigma_exp+',4)',),zb=zb,zs=zs,v_ice_rate=v_ice_rate)
phi = Expression('rho*g*(zs - x[0])*400*pow('+sigma_exp+',3)/(zs-zb)*slope',zb=zb,zs=zs,rho=rho,g=g,slope=slope)

# Load Interval mesh
mesh = IntervalMesh(n_mesh,zb,zs) 

# Create Solution Space
Q = FunctionSpace(mesh, "CG", 1)

# Define test and trial functions
theta, v = TrialFunction(Q), TestFunction(Q)

# Define last solution based on surface temperature
theta_last = interpolate(Expression("T0",T0=surface_temp(0)),Q)

# Variational Problem
a = (theta*v + dt*k/rho*cp*spy*dot(nabla_grad(theta),nabla_grad(v)) + spy*dt*dot(w,nabla_grad(theta))*v*dt)*dx
L = (theta_last*v)*dx + (qgeo*v*dt/rho/cp)*ds
#L = (theta_last*v - dt*u*dth_dx*v + phi*v/rho/cp*dt)*dx + (qgeo*v*dt/rho/cp)*ds
#L = (theta_last*v - dt*u*dth_dx*v + phi*v/rho/cp*dt)*dx + (qgeo*v*dt/rho/cp)*ds

# Boundary Conditions 
def DirichletBoundary(x, on_boundary):
  tol = 1E-14
  return on_boundary and (abs(x[0]-zs) < tol)

theta0 = Expression('-10 + 5*sin(2*pi*t/spy)',t=0,spy=spy)
bc = DirichletBC(Q, theta0, DirichletBoundary)

A = assemble(a)
b = None

theta = Function(Q)
t = dt
while t <= total_time:
  b = assemble(L, tensor=b)
  theta0.t = t
  bc.apply(A, b)
#  from IPython.core.debugger import Tracer
#  debug_here = Tracer()
#  debug_here()
  solve(A,theta.vector(), b)
  t += dt
  theta_last.assign(theta)

py.plot(theta.vector().array(),mesh.coordinates())
py.show()

