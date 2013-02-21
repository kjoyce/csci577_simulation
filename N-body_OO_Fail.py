#Object Oriented N-Body Problem
# I decided to jump in and flounder. Maybe you guys can see where my understanding of OO 
# begins and ends and we can discuss on friday. Hopefully this is in some way helpful to us all.


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
        return [x,v,m]

class System(object):
    def __init__(self,bodylist):
        self.n=len(bodylist)
    def state(self,bodylist):
        state=[]
        for i in range(n):
            state.append(bodylist[i].state) 
        state=array(state)
        return state                        
        
        
        
            
class Force(object):
    def __init__(self,x1,x2,m1,m2):
        self.x1=x1
        self.x2=x2
        self.m1=m1
        self.m2=m2
        G=6.67*10**(-11)
    def component_forces(self,x1,x2,m1,m2):
        r=norm(x2-x1)
        xdist=norm(x2[0]-x1[0])
        ydist=norm(x2[1]-x1[1])
        Fx=-G*m1*m2*xdist/r**3
        Fy=-G*m1*m2*ydist/r**3
        return array([Fx,Fy])
    
class Force_Pair(Force, System):
    def __init__(self):
        force_pairs=array([[0 for col in range(self.n)] for row in range(self.n)])
    
        for i in range(0,self.n):
            x1=System.state[i][0]
            m1=System.state[i][2]
            for j in range(i,self.n):
                x2=System.state[j][0]
                m2=System.state[j][2]
                force_pairs[i][j]=Force.component_forces(x1,x2,m1,m2)
                
