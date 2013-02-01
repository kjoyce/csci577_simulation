from scipy import *
from numpy import *
from coffee_data import time,black,cream

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
  
#Adjust parameter for black = .293, and cream = .279
def cooling_law(temp,c):
    return -c*(temp-23)

def eulers_method(f,y0,a,b,dx):    
    x = arange(a,b,dx)    
    y = ones(size(x))*y0
    for i in range(size(x)-1):
        y[i+1] = y[i] + f(x[i],y[i])*dx
    return (x,y)

#USAGE
# (t,mtblack) = eulers_method((lambda x,y: cooling_law(y,c_black)),82.3,0,46,.1)