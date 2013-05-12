"""
Simulation of Ice-Sheet temperature distribution
""" 
from dolfin import *
import numpy as np
import pylab as py
#  from IPython.core.debugger import Tracer
#  debug_here = Tracer()
#  debug_here()

spy = 31556926.

# Parameters with meaningful change 
# Search the comments for 'firn' to see what parameters to change 
# when doing firn experiment

# make contour plots
make_contour = False  # Change this to True when using firn analysis
contour_file_name = 'warming_contour.pdf'
ylims = (-2500,500)

zb = -2500. # Use this for realistic depth
#zb = -500. # Use this for seasonal variation firn
w_ice_rate = -.1
#w_ice_rate = 0 # Use this to turn off vertical advection
u_ice_rate = 30
#u_ice_rate = 0  # Use this to turn off horizontal advection
slope = .4E-5
#slope = 0  # Use this to turn off compression
dth_dx = 2E-5
qgeo = 60E-3
n_mesh = 40
#total_time = 60000*spy  # This should be long enough for two steady states
#total_time = 40*spy  #Use this for firn
total_time = 30000*spy  # Use this for one steady state
steady_state_tol = .0005
dt = 30*spy # Use this for steady state
#dt = .1*spy  # Use this for firn

# Constants 
zs = 10.
g = 9.81
rho = 600. # Use this for firn
rho = 911. # Use this for ice density
cp = 2009.
beta = 9.8E-8
theta_pmp = lambda z: -beta*rho*g*(zs - z)
k = 2.1
xlim = (-13,0)

# Surface temperature
def surface_temp(t):
  return -10 + 5*np.sin(2*pi*t/spy)

# Expressions
sigma_exp = '(x[0] - zb)/(zs - zb)'
sigma = Expression(sigma_exp,zb=zb,zs=zs)
w = Expression(('w_ice_rate*'+sigma_exp,),zb=zb,zs=zs,w_ice_rate=w_ice_rate)
u = Expression('u_ice_rate*pow('+sigma_exp+',4)',zb=zb,zs=zs,u_ice_rate=u_ice_rate)
phi = Expression('rho*g*(zs - x[0])*4*pow('+sigma_exp+',3)/(zs-zb)*slope',zb=zb,zs=zs,rho=rho,g=g,slope=slope)
bed_boundary = Expression('1-'+sigma_exp,zs=zs,zb=zb)

# Load Interval mesh
mesh = IntervalMesh(n_mesh,zb,zs) 

# Create Solution Space
Q = FunctionSpace(mesh, "CG", 1)

# Define test and trial functions
theta, v = TrialFunction(Q), TestFunction(Q)

# Define last solution based on surface temperature
theta_last = interpolate(Expression("T0",T0=surface_temp(0)),Q)
# Here is a bad hack that lets you start from a steady state
# theta_last.vector().set_local(array([ -2.19830112,  -2.48580721,  -3.03946921,  -3.58733924,
#        -4.12669701,  -4.65496409,  -5.16974199,  -5.66884429,
#        -6.1503219 ,  -6.61248099,  -7.0538933 ,  -7.47339887,
#        -7.87010168,  -8.24335848,  -8.592762  ,  -8.91811907,
#        -9.21942509,  -9.49683567,  -9.75063667,  -9.98121377,
#       -10.18902233, -10.37455858, -10.53833277, -10.6808447 ,
#       -10.80256222, -10.90390268, -10.98521761, -11.04678034,
#       -11.08877663, -11.11129785, -11.11433659, -11.0977842 ,
#       -11.06143005, -11.00496212, -10.92796855, -10.82994002,
#       -10.71027256, -10.56827068, -10.40315066, -10.21404386, -10.        ]))
# 
# Variational Problem
#a = (theta*v + dt*k/rho/cp*dot(nabla_grad(theta),nabla_grad(v)))*dx 
a = (theta*v + dt*k/rho/cp*dot(nabla_grad(theta),nabla_grad(v)) + dt*dot(w,nabla_grad(theta))*v/spy)*dx
#L = (theta_last*v - dt*u*dth_dx*v + phi*v/rho/cp*dt)*dx + (qgeo*v*dt/rho/cp)*ds
#L = (theta_last*v - dt*u*dth_dx*v + phi*v/rho/cp*dt)*dx + (qgeo*v*dt/rho/cp)*ds

# Boundary Conditions 
def DirichletBoundary(x, on_boundary):
  tol = 1E-14
  return on_boundary and (abs(x[0]-zs) < tol)

t_surface = Constant( surface_temp(0) )
bc = DirichletBC(Q, t_surface, DirichletBoundary)

A = assemble(a)
b = None

theta = Function(Q)
t = dt
py.ion()
py.figure(1)
py.clf()
ax, = py.plot(theta_last.vector().array(),mesh.coordinates(),'k.')
py.xlim(xlim)
py.xlabel("$\\theta$")
py.ylabel("$z$")

# for calculating steady state
steady_state_reached = False
last_ss_vector = py.ones(theta_last.vector().array().shape)*100 # so the loop runs make it different than before

# Data for contour plot
if make_contour:
  temp_contours = np.zeros((int(total_time/dt),len(mesh.coordinates())))
  i = 0

steady_state_time = {}  # key: time of steady state, value: temperature gradient
#while not(steady_state_reached):  # Use this for steady state experiment
while t <= total_time: # Use this for fixed time looping
### FEniCS Solving ###
  t_surface.assign( surface_temp(t) )
  #L = (theta_last*v)*dx + (bed_boundary*dt/rho/cp*qgeo*v)*ds
  #L = (theta_last*v - u*dt*dth_dx*v/spy)*dx + (bed_boundary*dt/rho/cp*qgeo*v)*ds
  L = (theta_last*v - dt*u*dth_dx*v/spy - dt*phi*v/rho/cp/spy)*dx + (bed_boundary*dt/rho/cp*qgeo*v)*ds

  b = assemble(L, tensor=b)
  bc.apply(A, b)
  solve(A,theta.vector(), b)
  t += dt
  theta_last = theta

#### Numpy calculations ####
  np_theta = theta_last.vector().array()
  # Check Bed Melting #
  z = mesh.coordinates()[:,0]
  idx = np_theta>theta_pmp(z)
  np_theta[idx] = theta_pmp(z)[idx]
  theta_last.vector().set_local(np_theta)

  # Check steady state by comparing every tenth temperature vector
  if (t%(10*dt) == 0):
    current_ss_vector = theta.vector().array()
    err = py.norm(current_ss_vector - last_ss_vector)/py.norm(last_ss_vector)
    steady_state_reached = err < steady_state_tol
    last_ss_vector = current_ss_vector
  if steady_state_reached:
    steady_state_time[str(t/spy)] = np_theta # key: time of steady state, value: temperature gradient
    def surface_temp(t): # Overload the old surface temperature
      return -8.
    steady_state_reached = False

  # Save values for contour map
  if make_contour:
    temp_contours[i,:] = np_theta
    i += 1
    
  py.title("Time: {}, Max temp: {:.3f}, Min Temp: {:.4f} PMP: {:.3f}".format(t/spy, max(np_theta), min(np_theta), theta_pmp(zb)))
  ax.set_xdata(np_theta)
  py.draw()


# Plot contour
def temp_contour_plot(filename,ylims):
  py.figure(figsize=(10,5))
  depth = mesh.coordinates()[:,0]
  py.ylim(ylims)
  py.contourf(np.arange(0,total_time/spy,dt/spy),depth,temp_contours.T,cmap=py.cm.winter) 
  py.colorbar()
  py.show()
  py.xlabel("Years") 
  py.ylabel("Elevation in Meters") 
  py.title("Temperature Over Time") 
  py.savefig(filename)
if make_contour:
  temp_contour_plot(contour_file_name,ylims)



