import matplotlib
import numpy
import math 
import scipy.optimize as sco
presence = numpy.array([[0,0,0],[0,3,0],[0,6,0],[0,9,0],[3,0,0],[3,3,0],[3,6,0],[3,9,0],[0,0,3],[0,3,3],[0,6,3],[0,9,3],[3,0,3],[3,3,3],[3,6,3],[3,9,3],
                        [0,0,6],[0,3,6],[0,6,6],[0,9,6],[3,0,6],[3,3,6],[3,6,6],[3,9,6],[0,0,9],[0,3,9],[0,6,9],[0,9,9],[3,0,9],[3,3,9],[3,6,9],[3,9,9],[0,0,0]])
target =   numpy.array([[0,0,9],[0,3,9],[0,6,9],[0,9,9],[3,0,9],[3,3,9],[3,6,9],[3,9,9],[0,0,6],[0,3,6],[0,6,6],[0,9,6],[3,0,6],[3,3,6],[3,6,6],[3,9,6],
                        [3,9,3],[3,6,3],[3,3,3],[3,0,3],[0,9,3],[0,6,3],[0,3,3],[0,0,3],[0,0,0],[3,9,0],[3,6,0],[3,3,0],[3,0,0],[0,9,0],[0,6,0],[0,3,0],[0,0,0]])
m = numpy.shape(presence)[0]
Kp = 10
Kc = 15
Dc = 1.5
Qm = 0
dt = 1e-2
r = numpy.zeros(m*m*3).reshape(m, m, 3)
p = numpy.zeros(m*m*3).reshape(m, m, 3)
for i in range(m):
    for j in range(m):
        r[i,j,:] = -target[i]+target[j]
for i in range(m):
    for j in range(m):
        p[i,j,:] = -presence[i]+presence[j]
temp = presence
tempv = numpy.zeros(m*3).reshape(m, 3)
log = presence
print(numpy.sum(numpy.abs(p-r)))
n = 0
def min_func_var(v):
    return numpy.dot(v,v)
while numpy.sum(numpy.sum(numpy.abs(p-r)))> 150:
    for i in range(m-1):    
        x0 = numpy.array([tempv[i][0],tempv[i][1],tempv[i][2],0])
        G = numpy.zeros((m+7)*4).reshape(m+7,4)
        [G[0][0],G[0][1],G[0][2]] = (r[i,m-1,:]-p[i,m-1,:])
        G[0][3] = -1
        for j in range(m-1):
            if i!= j:
                [G[j+1][0],G[j+1][1],G[j+1][2]]=  p[i,j,:]
        [G[m][0],G[m][1],G[m][2]]= numpy.array([1.0,0,0])       
        [G[m+1][0],G[m+1][1],G[m+1][2]]= numpy.array([0,1.0,0])
        [G[m+2][0],G[m+2][1],G[m+2][2]] =numpy.array([0,0,1])
        [G[m+3][0],G[m+3][1],G[m+3][2]]= numpy.array([-1.0,0,0])       
        [G[m+4][0],G[m+4][1],G[m+4][2]] = numpy.array([0,-1.0,0])
        [G[m+5][0],G[m+5][1],G[m+5][2]] = numpy.array([0,0,-1])
        G[m+6][3] = -1  
        h = numpy.zeros(m+7)
        a = (r[i,m-1,:]-p[i,m-1,:])
        h[0] = numpy.array([-(numpy.linalg.norm(a))**2*Kp])
        for j in range(m-1):
            if i!= j:
                h[j+1] = (Kc*(numpy.linalg.norm(p[i,j,:])**2-Dc**2)-Qm*numpy.linalg.norm(p[i,j,:]))
        h[m] = 5
        h[m+1] = 5
        h[m+2] = 5
        h[m+3] = 5
        h[m+4] = 5
        h[m+5] = 5
        h[m+6] = 0
        cons = [ {'type': 'ineq', 'fun': lambda x: h-numpy.dot(G, x)}]
        opts = sco.minimize(fun=min_func_var, x0=x0, method='SLSQP',constraints=cons)
        if opts['success'] == True:
            tempv[i] = numpy.array([opts['x'][0],opts['x'][1],opts['x'][2]])
        else:
            tempv[i] = tempv[i]*0.1
            print("check!!")
    tempv[m-1] = numpy.array([0,0,0])
    temp = temp + dt*tempv
    n = n+1
    if n%100 == 0:
        print(numpy.sum(numpy.abs(p-r)))
    for i in range(m):
        for j in range(m):
            p[i,j,:] = -temp[i]+temp[j]
    log = numpy.hstack((log,temp))
numpy.save("transitionlog6",log)