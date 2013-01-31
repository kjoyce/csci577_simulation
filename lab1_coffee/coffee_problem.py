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


