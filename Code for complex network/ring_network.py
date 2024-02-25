import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
import math
from tqdm import tqdm
import cmath
import csv


N = int(input("no. of nodes: "))
step = int(input("no. of time steps: "))
lmd_steps  = int(input("no. of lamda steps:"))

ring_matrix = [[ 0 for j in range(N)] for i in range(N)]
#print(ring_matrix)

for i in range(N-1):
    ring_matrix[i][i+1] = 1
    ring_matrix[i+1][i] = 1
ring_matrix[0][N-1] = 1
ring_matrix[N-1][0] = 1

ring_matrix = np.array(ring_matrix)
print(ring_matrix)


G = nx.from_numpy_matrix(ring_matrix)
#nx.draw(G,node_color='#8F00FF',with_labels = True, node_size=20)

#normal distribution
omega = np.zeros(N)
for i in range(N):
    omega[i] = (i)/(N-1)

    

alpha = 1


for node in G.nodes():
    G.nodes[node]['w'] = omega[node]
    G.nodes[node]['theta'] = np.random.rand()
    G.nodes[node]['thetadot'] = 0


nextG = G.copy() #copy graph

#nx.draw(G,node_color='#8F00FF', node_size=20) #to draw

#extracting adjacency matrix
p = nx.adjacency_matrix(G).toarray()
#print(p)

#changing weighted matrix to unweighted
A = (p>=1).astype(int)
print(A)
lmd0 = 0.02
final_order = []

"""Kuramoto Oscillator function"""
laus = []

for i in tqdm(range(lmd_steps)):
    lmd = lmd0*i
    laus.append(lmd)
    #print(laus)
    print(lmd)
    
    def f(t0,theta, node,t):
        cpl = 0


        # sigma terms
        for j in range(N):
            cpl += A[node][j]*np.sin(theta[j]-theta[node])*(abs(nextG.nodes[node]['w']-nextG.nodes[j]['w'])**alpha)


        thetadot = nextG.nodes[node]['w'] + (lmd/3)*cpl
        return thetadot

    # RK-4 method
    def rk4(t0,y0,xn,n):

        # Calculating step size
        h = (xn-t0)/n

        time_series = []
#        print(n,N)
        for i in range(n):
    #        data = np.zeros(N)
            temp = np.array(y0)

            time_series.append(y0) 
    #        print(time_series)
            for node in G.nodes():

                k1 = h * (f(t0, y0, node, h))
                k2 = h * (f((t0+h/2), (y0+k1/2),node, h))
                k3 = h * (f((t0+h/2), (y0+k2/2),node, h))
                k4 = h * (f((t0+h), (y0+k3),node, h))
                k = (k1+2*k2+2*k3+k4)/6

                yn = y0[node] + k

                temp[node] = yn
                t0 = t0+h
    #            data[node] = yn
            y0 = temp
        return time_series

    # Inputs
    t0 = 0
    y0 = np.zeros(N)
    for i in range(N):
        y0[i] = G.nodes[i]['theta']
    xn = 100
    
    

    data = np.array(rk4(t0,y0,xn,step))
    #data = np.array(data)

    def order(theta):
        z = sum(np.exp(theta*1j))/len(theta)
    #    print(np.absolute(z), np.angle(z))
        return np.absolute(z), np.angle(z)

    ordr = []
    angle = []

    st_ep = []
    for i in range(step):
        the_ta = data[i]
        x,y = order(the_ta)
        st_ep.append(i)
        ordr.append(x)
        angle.append(y)


    oreder = order(data[-1])
  
    final_order.append(oreder)
    print(oreder)

#print(ordr)
#plt.savefig("order.png")

#f.close()

plt.plot(laus,final_order[0:200])
plt.savefig("order2.png")