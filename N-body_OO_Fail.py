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
        return array([x],[v],m)

class System(object):
    def __init__(self,bodylist):
        self.n=len(bodylist)
    def state(self,bodylist):
        state=array([0 for col in range(2*n+n)])
        for i in range(n):
            state[2*i+1]=bodylist[i][0]     #position array x in first slot
            state[2*(i+1)]=bodylist[i][1]   #velocity array v in second slot
            state[-(n-i)]=bodylist[i][2]    #places masses at end of the list from m1,m2....mn
        return state                        #e.g. x1,v1,x2,v2,x3,v3,m1,m2,m3
                                            #Going deeper in: state=array([[x,y]1,[vx,vy]1,[x,y]2,
                                                            #[vx,vy]2,[x,y]3,[vx,vy]3,m1,m2,m3])
            
            
class Force(object):
    def __init__(self,x1,x2,m1,m2):
        self.x1=x1
        self.x2=x2
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
    #    x=System.state
        for i in range(0,self.n):
            x1=System.state[i][0]
            m1=System.state[i][2]
            for j in range(i,self.n):
                x2=System.state[j][0]
                m2=System.state[j][2]
                force_pairs[i][j]=Force.component_forces(x1,x2,m1,m2)
                
