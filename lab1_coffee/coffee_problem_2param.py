from scipy import *
from numpy import *
from coffee_data import time,black,cream
from IPython.core.debugger import Tracer
debug_here = Tracer()

# Transform the data 
def transform_data(temp,T0):
  return log(temp - 23)

tblack = transform_data(black,black[0])
tcream = transform_data(cream,cream[0])

# Solve the linear least squares problem
# A c = temp + eps
# is (A'A)^(-1) A' temp
def lst_square(time,data):
  n = size(time)
  A = column_stack((ones(n), time))
  return dot( linalg.inv( dot(A.T, A)) , dot( A.T, data) )

# the minus is so c > 0
black_param = lst_square( -time, tblack )
cream_param = lst_square( -time, tcream )

def actual_T(t,c,T0):
  return (T0 - 23)*exp(-c*t) + 23


