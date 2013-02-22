#Object Oriented N-Body Problem


#Body object to create n bodies with, which will be passed into System via bodylist
class Body(object):
    x=[0,0]
    v=[0,0]
    m=0
    def __init__(self,x,v,m):
        self.x=x
        self.v=v
        self.m=m
        
    def state(self):
        return [self.x,self.v,self.m]


class System(object):
    def __init__(self,bodylist):
        self.bodylist=bodylist
        self.n=len(bodylist)
    def get_state(self):
        sys_state=[]
        for i in range(self.n):
            sys_state.append(self.bodylist[i].state())
        
        return sys_state
        
        
        
            
class Force(System):
    def __init__(self,system):
        self.n=system.n
        self.system=system
        self.G=6.67*10**(-11)
    def component_forces(self,x1,x2,m1,m2):
        r=norm(x2-x1)
        xdist=norm(x2[0]-x1[0])
        ydist=norm(x2[1]-x1[1])
        Fx=-self.G*m1*m2*xdist/r**3
        Fy=-self.G*m1*m2*ydist/r**3
        return array([Fx,Fy])
    
    def net_force(self):
        
        
        force_pairs=array([[0 for col in range(self.n)] for row in range(self.n)])
        Fxy_net=[]
        
        for i in range(self.n):
            x1=self.system.get_state()[i][0]
            m1=self.system.get_state()[i][2]
            for j in range(i+1,self.n):
                x2=self.system.get_state()[j][0]
                m2=self.system.get_state()[j][2]
                force_pairs[i+1][j]=self.component_forces(self.x1,self.x2,self.m1,self.m2)
        
        for i in range(self.n):
            sum_force=sum(force_pairs, axis=1)[i]+sum(force_pairs, axis=0)[i]-force_pairs[i][i+1]
            Fxy_net.append(sum_force)
        
        return Fxy_net
        

