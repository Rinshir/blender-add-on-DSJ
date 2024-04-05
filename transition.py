import matplotlib
import numpy
import matplotlib.pyplot as plt
import math 
import cvxopt
from cvxopt import matrix
presence = numpy.array([[0,3],[0,0],[3,3],[3,0],[6,3],[6,0],[9,3],[9,0]])
target = numpy.array([[0,9],[0,6],[0,3],[0,0],[3,9],[3,6],[3,3],[3,0]])
Kp = 100
Kc = 1
Dc = 0.1
Qm = 0.01
dt = 1.0
r = numpy.zeros(128).reshape(8, 8, 2)
p = numpy.zeros(128).reshape(8, 8, 2)
for i in range(8):
    for j in range(8):
        r[i,j,:] = target[i]-target[j]
for i in range(8):
    for j in range(8):
        p[i,j,:] = presence[i]-presence[j]
v = numpy.zeros(16).reshape(8, 2,1)
temp = presence
tempv = numpy.zeros(16).reshape(8, 2)
log = presence
print(p-r)
while numpy.sum(numpy.abs(p-r)) > 1:
    for i in range(8):
        A=matrix(numpy.array([[1.0,0.0],[0.0,1.0]]))
        q=matrix(numpy.array([0.0,0.0]))
        G = numpy.zeros(18).reshape(9,2)
        h = numpy.zeros(9)
        G[0] = numpy.sum(r[i,:,:]-p[i,:,:],axis=0)
        h = numpy.zeros(9)
        for j in range(8):
            if i!= j:
                    G[j+1]= p[i,j,:]
        G = matrix(G)
        for j in range(8):
            h[0] = numpy.array([-numpy.linalg.norm(numpy.sum(r[i,:,:]-p[i,:,:],axis=0))**2*Kp])
            if i!= j:
                h[j+1] = (Kc*(numpy.linalg.norm(p[i,j,:])**2-Dc**2)-Qm*numpy.linalg.norm(p[i,j,:]))
        h = matrix(h)
        sol=cvxopt.solvers.qp(A,q,G,h)
        tempv[i] = numpy.array([sol["x"][0],sol["x"][1]])
    for i in range(8):
        for j in range(8):
            p[i,j,:] = temp[i]-temp[j]
    temp = temp + dt*tempv
    log = numpy.hstack((log,temp))
    plt.scatter(temp[:,0],temp[:,1])
    print(p-r)
plt.scatter(temp[:,0],temp[:,1])
plt.savefig("temp.png")