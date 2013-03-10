Created on Fri Mar  8 22:43:52 2013

@author: will

from pylab import *

from matplotlib import animation
from IPython.core.debugger import Tracer
debug_here=Tracer()
debug_here()


epsilon=1
sigma=1
mass=1
    

class container(object):
    def __init__(self,dimensions):
        self.Lx=dimensions[0]
        self.Ly=dimensions[1]
        self.Lz=dimensions[2]
        self.x=array([])
        self.vx=array([])
        self.y=array([])
        self.vy=array([])
        self.z=array([])
        self.vz=array([])
        self.mass=array([])
        self.ax,self.ay,self.az =self.acceleration()
        
    
    def addParticle(self,x,y,z,vx,vy,vz,m):
        self.x=hstack((self.x, x))
        self.y=hstack((self.y, y))      
        self.z=hstack((self.z, z))
        self.vx=hstack((self.vx, vx))
        self.vy=hstack((self.vy, vy))
        self.vz=hstack((self.vz, vz))
        self.mass=hstack((self.mass, m))
        self.ax,self.ay,self.az=self.acceleration()
        
        
    def d_matrix(self,A,L):
        
        A=array(A)
        Asquare=tile(A, [size(A),1])
        if L == 0:
            Adist=(Asquare-Asquare.T)
        else:
            Adist=(Asquare-Asquare.T)%L
            Adist[Adist>L/2.]=Adist[Adist>L/2.]-L
            Adist[Adist<-L/2.]=Adist[Adist<-L/2.]+L
        
        return Adist
        
    def force_I(self,component,r):
        mag=(24.*epsilon/r)*(2*(sigma/r)**12-(sigma/r)**6)
        f_matrix=(mag*component/r)
        f_matrix=nan_to_num(f_matrix)        
        
        f_net=sum(f_matrix,axis=1)
        return f_net
        
    
    def acceleration(self):
        xdist=self.d_matrix(self.x, self.Lx)
        ydist=self.d_matrix(self.y, self.Ly)
        zdist=self.d_matrix(self.z, self.Lz)
        
        r=sqrt(xdist**2+ydist**2+zdist**2)        
        
        ax=self.force_I(xdist,r)/self.mass
        ay=self.force_I(ydist,r)/self.mass
        az=self.force_I(zdist,r)/self.mass
        self.ax,self.ay,self.az=ax,ay,az
        return self.ax,self.ay,self.az

class integrator(object):
    def __init__(self,c,dt):
        self.c=c
        self.dt=dt
        
        
    def verlet(self):
        c=self.c
        dt=self.dt
        
        old_ax,old_ay,old_az=c.ax,c.ay,c.az
        c.x=c.x+c.vx*dt+.5*c.ax*dt**2
        c.y=c.y+c.vy*dt+.5*c.ay*dt**2
        c.z=c.z+c.vz*dt+.5*c.az*dt**2
        c.ax,c.ay,c.az=c.acceleration()
        c.vx=c.vx+.5*(c.ax+old_ax)*dt
        c.vy=c.vy+.5*(c.ay+old_ay)*dt
        c.vz=c.vz+.5*(c.az+old_az)*dt



lx=10.
ly=10.
lz=0
dims=[lx,ly,lz]
n=8.
init_vel=.5
masses=1.
radius=2.**(1./6.)

xpos=array(linspace(radius,lx-radius,n))
ypos=zeros(n)+ly/2
zpos=zeros(n)
xvel=zeros(n)
yvel=zeros(n)+init_vel
zvel=zeros(n)

#initialize container for line test

class containerinit(object):
    def __init__(self,container,initialization):
        self.c=container
        c=self.c
        c.Lx=10.
        dist = c.Lx / 5.
        vel = dist /5.
        
        if initialization == 'one':
            c.addParticle(0,dist,0,0,0,0,1)
 
        elif initialization == 'two':
            c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
            c.addParticle(dist,0.,0,-vel,0.,0.,1.)
 
        elif initialization == 'three':
            c.addParticle(0.,dist*sqrt(3)/2,0.,0.,-vel,0.,1.)
            c.addParticle(-dist,0.,0.,vel,0.,0.,1.0)
            c.addParticle(dist,0.,0,-vel,0.,0.,1.)
 
        elif initialization == 'four':
            c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
            c.addParticle(dist,0.,0,-vel,0.,0.,1.)
            c.addParticle(0.,dist,0,0.,-vel,0.,1.)
            c.addParticle(0.,-dist,0,0.,vel,0.,1.)
            
        elif initialization == 'six':
            c.addParticle(0.,dist,0.,0.,-vel,0.,1.)
            c.addParticle(0.,-dist,0.,0.,vel,0.,1.)
            c.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
            c.addParticle(-dist/sqrt(2),dist/sqrt(2),0,vel/sqrt(2),-vel/sqrt(2),0.,1.)
            c.addParticle(-dist/sqrt(2),-dist/sqrt(2),0,vel/sqrt(2),vel/sqrt(2),0.,1.)
            c.addParticle(dist/sqrt(2),-dist/sqrt(2),0,-vel/sqrt(2),vel/sqrt(2),0.,1.)
            
        elif initialization == 'eight':
            c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
            c.addParticle(dist,0.,0,-vel,0.,0.,1.)
            c.addParticle(0.,dist,0,0.,-vel,0.,1.)
            c.addParticle(0.,-dist,0,0.,vel,0.,1.)
 
            c.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
            c.addParticle(-dist/sqrt(2),dist/sqrt(2),0,vel/sqrt(2),-vel/sqrt(2),0.,1.)
            c.addParticle(-dist/sqrt(2),-dist/sqrt(2),0,vel/sqrt(2),vel/sqrt(2),0.,1.)
            c.addParticle(dist/sqrt(2),-dist/sqrt(2),0,-vel/sqrt(2),vel/sqrt(2),0.,1.)
        
        elif initialization == 'line':
            
            N=8.
            init_vel=.5
            masses=1.
            radius=2.**(1./6.)
            xpos=array(linspace(radius,c.Lx-radius,N))
            ypos=c.Ly/2
            zpos=0
            xvel=0
            yvel=init_vel
            zvel=0
            for i in range(int(n)):
                c.addParticle(xpos[i],ypos[i],zpos,xvel,yvel,zvel,masses)
        
        elif initialization == 'square_lattice':
            N = 8             # Particles per row
            c.Ly = c.Lx       # Extents determined by Lx input
            d = 2.**(1/6.)    # Particle diameter
            x = linspace(-c.Lx/2+d/2.,c.Lx/2-d/2,N)
            y = linspace(-c.Lx/2+d/2.,c.Lx/2-d/2,N)
            for i in range(x.size):
                for j in range(y.size):
                    c.addParticle(x[i],y[j],0,0,0,0,1)
        
#         elif initialization == 'triangular_lattice':
#             N = 8                       # particles per row
#             c.Ly = sqrt(3) / 2. * c.Lx  # Set this based on Lx
#             d = 2.**(1/6.)              # diameter
#             x =  linspace(-c.Lx/2 + 3.*d/4.,c.Lx/2. - 1.*d/4., N) # Unstaggered
#             xs = linspace(-c.Lx/2 + d/4.   ,c.Lx/2. - 3.*d/4., N) # Staggered
#             y =  linspace(-c.Ly/2 + d/2.,c.Ly/2  - d/2, N)
# 
#             for i in range(N):
#                 for j in range(N):
#                     if mod(i,2)==0:
#                         c.addParticle(x[j],y[i],0,0,0,0,1)
#                     else:
#                         c.addParticle(xs[j],y[i],0,0,0,0,1)   



lx=10.
ly=10.
lz=0
dims=[lx,ly,lz]
dt=.01
c=container(dims)

  
c=container(dims)
containerinit(c,'four')

    
    
    
    


#animation

forward=integrator(c,dt)

circles=[]
fig=plt.figure(figsize=(6,6))
ax=plt.gca()
ax.set_aspect('equal')
ax.set_xlim((0,c.Lx))
ax.set_ylim((0,c.Ly))

def prettify_circle(e):
    
    color="lightsteelblue"
    facecolor="green"
    alpha=.6
    e.set_clip_box(ax.bbox)
    e.set_edgecolor( color )
    e.set_linewidth(3)
    e.set_facecolor( facecolor )  # "none" not None
    e.set_alpha( alpha )
    return e
  
for i in range(0,size(c.x)):
    e = Circle( (c.x[i],c.y[i]), radius=.5*2.**(1./6.), animated=True)
    e = prettify_circle(e)
    circles.append(ax.add_patch(e))

def init():
    
    return circles

def next_frame(i):
#  print "Frame: {}".format(i)
    forward.verlet()
    for i in range(len(circles)):
        e.center=(c.x[i],c.y[i])
    return circles


anim = animation.FuncAnimation(fig,next_frame,init_func=init,\
frames=200, interval=1, blit=True)

plt.show()


