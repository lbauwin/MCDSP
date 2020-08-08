# # Minimum connected dominating set Problem

import numpy as np
import matplotlib.pyplot as plt
from docplex.mp.model import Model
from MDS import *
from SSL import *
from SSL_lazy import *
from Martin import *
from MTZ import *
from SCF import *
import sys

rnd = np.random
#rnd.seed(3)


# Fuctions to define graphs
def random_graph(v,e):
    #First edges to assure graph is connected
    initialSet = []
    visitedSet = []
    edges =[]
    vertices = []
    for i in range(v):
        initialSet.append(i+1)
        vertices.append(i+1)
    curVertex = rnd.choice(initialSet)
    initialSet.remove(curVertex)
    visitedSet.append(curVertex)
    edgeCnt = 0
    while initialSet:
        adjVertex = rnd.choice(initialSet)
        edge = (curVertex, adjVertex)
        edges.append(edge)
        edgeCnt+=1
        initialSet.remove(adjVertex)
        visitedSet.append(adjVertex)
        curVertex=adjVertex
    #Second add all other edges randomly
    while edgeCnt < e:
        node1 = rnd.randint(1,v+1)
        node2 = rnd.randint(1,v+1)
        iteration = 0
        while (node1 == node2 or (node1,node2) in edges or (node2,node1) in edges): #If same random number
            node2 += 1
            if (node2 == v+1 ):
                node2 = 1
            if iteration == v:
                node1 +=1
                if (node1==v+1):
                    node1 = 1
                iteration = 0
            iteration+=1
        edges.append((node1,node2))
        edgeCnt+=1
    return edges


def IEEE_14_Bus_graph():
    return [(1,2),(1,5),(2,3),(2,4),(2,5),(3,4),(4,5),(4,7),(4,9),(5,6),(6,11),(6,12),(6,13),(7,8),(7,9),
         (9,10),(9,14),(10,11),(12,13),(13,14)]

def IEEE_30_Bus_graph():
    return [(1,2),(1,3),(2,4),(2,5),(2,6),(3,4),(4,6),(4,12),(5,7),(6,7),(6,8),(6,9),(6,10),(6,28),(8,28),(9,10),
        (9,11),(10,17),(10,20),(10,21),(10,22),(12,13),(12,14),(12,15),(12,16),(14,15),(15,18),(15,23),(16,17),
        (18,19),(19,20),(21,22),(22,24),(23,24),(24,25),(25,26),(25,27),(27,28),(27,29),(27,30),(29,30)]

def IEEE_57_Bus_graph():
    return [(1,2),(1,16),(1,17),(1,15),(2,3),(3,4),(3,15),(4,5),(4,6),(4,18),(5,6),(6,7),(6,8),(7,8),(7,29),(8,9),(9,10),(9,55),(9,11),(9,13),(9,12),(10,12),(10,51),(11,13),(11,41),(11,43),
        (12,13),(12,16),(12,17),(13,14),(13,15),(13,49),(14,15),(14,46),(15,45),(18,19),(19,20),(20,21),(21,22),(22,23),(22,38),(23,24),(24,25),(24,26),(25,30),(26,27),(27,28),
        (28,29),(29,52),(30,31),(31,32),(32,33),(32,34),(34,35),(35,36),(36,37),(36,40),(37,38),(37,39),(38,44),(38,48),(38,49),(39,57),(40,56),(41,42),(41,43),(41,56),(42,56),(44,45),
        (46,47),(47,48),(48,49),(49,50),(50,51),(52,53),(53,54),(54,55),(56,57)]



def main():
    #get argument
    if(len(sys.argv)<2):
        print("Use the command as:\npython3 Tests.py <method> <nbr of vertices> <degree of nodes>\nThe method can be:\n0 \trandomized graph\n1 \tIEEE_14_Bus\n2 \tIEEE_30_Bus\n3 \tIEEE_57_Bus")
        sys.exit(1)
    elif (sys.argv[1]=="0" and len(sys.argv)<4):
        print("Need the number of vertices and edges for a random graph.\nUse the command as:\npython3 Tests.py <method> <nbr of vertices> <degree of nodes>")
        sys.exit(1)
    else:
        switcher={
                0:'random',
                1:'IEEE-14-Bus',
                2:'IEEE-30-Bus',
                3:'IEEE-57-Bus',
             }
        name= switcher.get(int(sys.argv[1]), 0)
        if (name==0):
            print("The method available are:\n0 \trandomized graph\n1 \tIEEE_14_Bus\n2 \tIEEE_30_Bus\n3 \tIEEE_57_Bus")
            sys.exit(1)

    E = [] # Set of vertices
    if name=="random":
        v = int(sys.argv[2]) # Number of nodes
        e = int(float(sys.argv[3])*v/2) # Formula to get average degree
        E = random_graph(v,e)

    elif name == "IEEE-14-Bus":
        v=14
        E = IEEE_14_Bus_graph()
    elif name == "IEEE-30-Bus":
        v=30
        E = IEEE_30_Bus_graph()

    elif name == "IEEE-57-Bus":
        v=57
        E = IEEE_57_Bus_graph()

    V = [i for i in range(1,v+1)] # Set of vertices
    A = [] # Neighboorhood matrix
    for i in range(v):
        A.append([])
        for j in range(v):
            if (i+1,j+1) in E or (j+1,i+1) in E or i==j:
                A[i].append(1)
            else:
                A[i].append(0)

    #orig_stdout = sys.stdout
    #filename = 'out_'+name+"_"+str(v)+"_"+str(len(E))+'.txt'
    #f = open(filename, 'a')
    #sys.stdout = f
    # start solver
    status = True
    #MDS(V,E,A, status)
    print("\n\nSolving MTZ...")
    Miller_Tucker_Zemlin(V,E,A, status)
    print("\n\nSolving SSL...")
    Simonetti_SallesDaCunha_Lucena(V,E,A, status)
    Simonetti_SallesDaCunha_Lucena_Lazy(V,E,A, status)

    print("\n\nSolving SCF...")
    Single_Commodity_Flow(V,E,A, status)
    print("\n\nSolving Martin...")
    Martin(V,E,A, status)
    #print("\n\n----------------------------------------------------\n\n")
    #sys.stdout = orig_stdout
    #f.close()


if __name__ == "__main__":
    main()
