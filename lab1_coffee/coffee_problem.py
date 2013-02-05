from scipy import *
from numpy import *
from coffee_data import time,black,cream
from matplotlib.pyplot import *
from matplotlib.mlab import find

# Transform the data 
def transform_data(temp,T0):
  return log((temp - 23)/(T0 - 23))

tblack = transform_data(black,black[0])
tcream = transform_data(cream,cream[0])

# Solve the linear least squares problem
# A c = temp + eps
# is (A'A)^(-1) A' temp
c_black = 1/dot(time.T,time) * dot(-time.T,tblack)
c_cream = 1/dot(time.T,time) * dot(-time.T,tcream)

def actual_T(t,c,T0):
  return (T0 - 23)*exp(-c*t) + 23
  
def cooling_law(temp,c):
    return -c*(temp-23)

def eulers_method(f,y0,a,b,dx):    
    x = arange(a,b,dx)    
    y = ones(size(x))*y0
    for i in range(size(x)-1):
        y[i+1] = y[i] + f(x[i],y[i])*dx
    return (x,y)

# Evaluate convergence of Euler's Method for black coffee
(t,mtblack) = eulers_method((lambda x,y: cooling_law(y,c_black)),black[0],0,46,.1)
mblack = actual_T(t,c_black,black[0])
err = max(abs(mtblack - mblack))/max(abs(mblack)) 
print '% error for Euler\'s method: {0:4f}'.format(err)

# use Euler's method for models
(t,mtblack) = eulers_method((lambda x,y: cooling_law(y,c_black)),black[0],0,46,.1)
(t,mtcream) = eulers_method((lambda x,y: cooling_law(y,c_cream)),cream[0],0,46,.1)

figure(1)
clf()
plot(t,mtblack,'-b')
plot(time,black,'.b')
plot(t,mtcream,':r')
plot(time,cream,'xr')
legend(['model fit black','black data','model fit cream','cream data']) 
title('Model fits')
xlabel('time (minutes)')
ylabel('temperature ($C^0$)')
show()


# use Euler's method to answer question
(t,test_black) = eulers_method((lambda x,y: cooling_law(y,c_black)),90,0,46,.1)
(t,test_cream) = eulers_method((lambda x,y: cooling_law(y,c_cream)),85,0,46,.1)
last_idx = find(test_black <= 80)[0]
last_idxx = find(test_cream <= 75)[0]
# Compare t[last_idx] to t[last_idxx]

figure(2)
clf()
ax = gca()
plot(t[0:(last_idx+1)],test_black[0:(last_idx+1)])

plot(t[0:(last_idxx+1)],test_cream[0:(last_idxx+1)],'r')

ylim((70,95))
l = Line2D([t[last_idx],t[last_idx]],[0,test_black[last_idx]])
l.set_linestyle(':') 
l.set_color('blue')
ax.add_line(l)
hl = Line2D([0,t[last_idx]],[test_black[last_idx],test_black[last_idx]]);
hl.set_linestyle(':')
hl.set_color('blue')
ax.add_line(hl)


ll = Line2D([t[last_idxx],t[last_idxx]],[0,test_cream[last_idxx]])
ll.set_linestyle(':') 
ll.set_color('red')
ax.add_line(ll)
hll = Line2D([0,t[last_idxx]],[test_cream[last_idxx],test_cream[last_idxx]]);
hll.set_linestyle(':')
hll.set_color('red')
ax.add_line(hll)
legend(['black coffee temp.','creamed coffee temp.']) 
title('Temperature Comparison')
xlabel('time (minutes)')
ylabel('temperature ($C^0$)')
show()

