class Container(object):
    def __init__(self,x,vx,xdim,y,vy,ydim,z,vz,zdim):
        self.x,self.vx,self.xdim=x,vx,xdim
        self.y,self.vy,self.ydim=y,vy,ydim
        self.z,self.vz,self.zdim=z,vz,zdim
    
    def d_matrix(self):
        x=tile(self.x, size(self.x),1)
        y=tile(self.y, size(self.y),1)
        z=tile(self.z, size(self.z),1)
        xdim=self.xdim
        ydim=self.ydim
        zdim=self.zdim
        xdist=(x-x.T)%xdim
        ydist=(y-y.T)%ydim
        zdist=(z-z.T)%zdim
        xdist[xdist>xdim/2.]=xdist[xdist>xdim/2.]-xdim
        xdist[xdist<-xdim/2.]=xdist[xdist<-xdim/2.]+xdim
        ydist[ydist>ydim/2.]=ydist[ydist>ydim/2.]-ydim
        ydist[ydist<-ydim/2.]=ydist[ydist<-ydim/2.]+ydim
        zdist[zdist>zdim/2.]=zdist[zdist>zdim/2.]-zdim
        zdist[zdist<-zdim/2.]=zdist[zdist<-zdim/2.]+zdim
        r=sqrt(xdist**2+ydist**2+zdist**2)
        self.xdist,self.ydist,self.zdist,self.r=xdist,ydist,zdist,r
        return xdist, ydist, zdist, r
        
    def force_I(self,component,r):
        f_matrix=24.*(2*(1./r**12)-(1./r)**6)*component/r**2
        f_net=sum(f_matrix,axis=1)
        return f_net
        
    
    def __call__(self):
        xdist, ydist, zdist, r = d_matrix()
        ax=force_I(xdist,r)
        ay=force_I(ydist,r)
        az=force_I(zdist,r)
        self.ax,self.ay,self.az=ax,ay,az
        return self.x,self.y,self.z,self.ax,self.ay,self.az

class integrator(object):
    def __init__(self,c,dt):
        self.c=c
        self.__dt=dt
        
        
    def int_I(self):
        c=self.c
        dt=self.__dt
        c.x=c.x+c.vx*dt+.5*c.ax*dt**2
        c.y=c.y+c.vy*dt+.5*c.ay*dt**2
        c.z=c.z+c.vz*dt+.5*c.az*dt**2
        new=Container(c.x,c.vx,c.xdim,c.y,c.vy,c.ydim,c.z,c.vz,c.zdim)
        c.vx=c.vx+.5*(c.ax+new.ax)*dt
        c.vy=c.vy+.5*(c.ay+new.ay)*dt
        c.vz=c.vz+.5*(c.az+new.az)*dt
