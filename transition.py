import matplotlib
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math 
import cvxopt
from cvxopt import matrix, solvers
presence = numpy.array([[0,3],[0,0],[3,3],[3,0],[6,3],[6,0],[9,3],[9,0]])
target = numpy.array([[0,9],[0,6],[0,3],[0,0],[3,9],[3,6],[3,3],[3,0]])
Kp = 1000
Kc = 100000
Dc = 0.1
Qm = 0
dt = 1.0e-8
r = numpy.zeros(128).reshape(8, 8, 2)
p = numpy.zeros(128).reshape(8, 8, 2)
for i in range(8):
    for j in range(8):
        r[i,j,:] = -target[i]+target[j]
for i in range(8):
    for j in range(8):
        p[i,j,:] = -presence[i]+presence[j]
v = numpy.zeros(16).reshape(8, 2,1)
temp = presence
tempv = numpy.zeros(16).reshape(8, 2)
log = presence
print(numpy.sum(numpy.abs(p-r)))
fig  = plt.figure
ims = []
while numpy.sum(numpy.abs(p-r)) > 10:
    for i in range(8):
        A=matrix(numpy.array([[1.0,0.0],[0.0,1.0]]))
        q=matrix(numpy.array([0.0,0.0]))
        G = numpy.zeros(18).reshape(9,2)
        h = numpy.zeros(9)
        for j in range(8):
            G[0] = G[0]+(r[j,i,:]-p[j,i,:])
        h = numpy.zeros(9)
        for j in range(8):
            if i!= j:
                    G[j+1]= p[j,i,:]
        G = matrix(G)
        for j in range(8):
            for k in range(8):
                a = 0
                a = a+(r[i,k,:]-p[i,k,:])
            h[0] = numpy.array([-(numpy.linalg.norm(a))**2*Kp])
            if i!= j:
                h[j+1] = (Kc*(numpy.linalg.norm(p[i,j,:])**2-Dc**2)-Qm*numpy.linalg.norm(p[i,j,:]))
        h = matrix(h)
        tempv[i,0]= solvers.coneqp(A,q,G,h)["x"][0]
        tempv[i,1]= solvers.coneqp(A,q,G,h)["x"][1]
    for i in range(8):
        for j in range(8):
            p[i,j,:] = temp[i]-temp[j]
    temp = temp + dt*tempv
    log = numpy.hstack((log,temp))
    im = plt.scatter(temp[:,0],temp[:,1])
    ims.append(im)
    plt.gca().clear()
ani = animation.ArtistAnimation(fig,ims,interval= 100)
ani.save("transition.mp4")