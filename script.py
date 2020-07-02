import sys
import os

def benchmark(iteration, algorithm, nodes=None, edges=None):
    for i in range(iteration):
        if algorithm==0:
            os.system("python3 Tests.py "+str(algorithm)+" "+str(nodes)+" "+str(edges))
        else:
            os.system("python3 Tests.py "+str(algorithm))

if __name__=='__main__':
    #benchmark(10, 1) done
    #benchmark(10, 2) done
    #benchmark(3, 3)
    nodes = [200]
    density = [70,60,50,40]
    for i in range(len(nodes)):
        for j in range(len(density)):
            E = nodes[i]*(nodes[i]-1)*density[j]/200
            edges = int(round(2*E/nodes[i]))
            benchmark(10, 0, nodes[i], edges)
