import pickle
import re
from numpy import array,ones,dot,arange
from numpy.linalg import solve
import matplotlib.pyplot as plt
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
#debug_here()

def lst_square_params(xy):
  A = ones(xy.shape)
  A[:,1] = xy[:,0]
  y = xy[:,1]
  return solve(dot(A.T,A),dot(A.T,y))

def slicedict(d, s): 
  return {int(re.search('(?<=load)-?\d+',k).group(0)):v for k,v in d.iteritems() if k.startswith(s)} 

#f = open("2013-03-24_21:18force_run.dump","r")
f = open("2013-03-25_01:27force_run.dump","r") 
force,out = pickle.load(f)

sld9  = array(slicedict(force,'sled9').items())
yhat9 = lst_square_params(sld9)
sld13 = array(slicedict(force,'sled13').items())
yhat13 = lst_square_params(sld13)
sld17 = array(slicedict(force,'sled17').items())
yhat17 = lst_square_params(sld17)

plt.figure()
x = arange(-20,40)
plt.plot(x,yhat9[0] + x*yhat9[1],'b:')
plt.plot(x,yhat13[0] + x*yhat13[1],'r:')
plt.plot(x,yhat17[0] + x*yhat17[1],'g:')
plt.plot(sld9[:,0],sld9[:,1],'bo')
plt.plot(sld13[:,0],sld13[:,1],'r*')
plt.plot(sld17[:,0],sld17[:,1],'g+')
plt.xlabel("$W$ - Average Load per particle (non-dimensionalized)")
plt.ylabel("Maximum Pull Force (non-dimensionalized)")
plt.title("Average Load vs. Max Pull Force")
plt.legend(plt.gca().lines[-3:],("9-Configuration","13-Configuration","17-Configuration"),loc="best")
plt.savefig("load_maxforce.pdf")
plt.show()

print "beta0: {} beta1: {}".format(*yhat9)
print "beta0: {} beta1: {}".format(*yhat13)
print "beta0: {} beta1: {}".format(*yhat17)
