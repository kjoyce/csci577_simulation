from numpy import array,size,tile
#from IPython.core.debugger import Tracer

class DistMatrix(object):
    def __init__(self):
        pass

    def CalcDist(self, x, y=[], z=[]):
        tilex = tile(x , (size(x),1))
        distx = tilex.T - tilex
        tiley = tile(y , (size(y),1))
        disty = tiley.T - tiley
        tilez = tile(z , (size(z),1))
        distz = tilez.T - tilez
        return distx, disty, distz

    def UpdatePBoundries(self, vector, domain):
        xupdateBoundries = vector
        xupdateBoundries % domain
        return xupdateBoundries

    def CheckDist(self, distMatrix, domain):
        checkDist = distMatrix
        checkDist[checkDist > domain/2.] %= domain
        checkDist[checkDist < -domain/2.] %= domain
        return checkDist

matrix = DistMatrix()
x1 = array([.5,2.5,4.,2.])
y1 = array([4,5,8,4])
matrix.UpdatePBoundries(x1,3.)
matrix.UpdatePBoundries(y1,3.)
test1,test2, test3 = matrix.CalcDist(x1,y1)

print test1
matrix.CheckDist(test1,3.)

print x1
print test1
print test2
