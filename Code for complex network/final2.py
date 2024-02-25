

import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
import math
from tqdm import tqdm
import cmath



N =100

G = nx.scale_free_graph(N) #scale free directed network is created
G = G.to_undirected() # changed to undirected graphG.pos = nx.spring_layout(G) #to give nodes spring like structure
    #nx.draw(G,node_color='#8F00FF', node_size=20) #to draw

#normal distribution
omega = np.random.normal(loc=1, scale=2, size=(N))


for node in G.nodes():
    G.nodes[node]['w'] = omega[node]
    G.nodes[node]['theta'] = np.random.rand()
    G.nodes[node]['thetadot'] = 0

nextG = G.copy() #copy graph

nx.draw(G,node_color='#8F00FF', node_size=20) #to draw

#extracting adjacency matrix
p = nx.adjacency_matrix(G).toarray()
#print(p)
#changing weighted matrix to unweighted
A = (p>=1).astype(int)
print(A)
lmd0 = 0.01
final_order = []


"""Kuramoto Oscillator function"""
laus = []
for i in tqdm(range(500)):
    lmd = 0.5 + lmd0*i
    laus.append(lmd)
    print(laus)
    print(lmd)
    
    def f(t0,theta, node,t):
        cpl = 0


        # sigma terms
        for j in range(N):
            cpl += A[node][j]*np.sin(theta[j]-theta[node])


        thetadot = nextG.nodes[node]['w'] + (lmd*abs(nextG.nodes[node]['w'])*cpl)/(sum(A[node]))
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
    step = 20000 

    data = np.array(rk4(t0,y0,xn,step))
    #data = np.array(data)
#    print(np.shape(data))

#    print(data)

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

    #ordr = np.arra

    oreder = order(data[-1])
    final_order.append(oreder)
    print(oreder)
#    print(ordr)
    #print(ordr)
    #print(angle)
    #plt.plot(st_ep,ordr)
    #plt.plot(angle)
    #plt.savefig("order.png")
plt.plot(laus,final_order[0:200])
plt.savefig("order50.png")
