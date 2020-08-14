from docplex.mp.model import Model
import numpy as np
from time import time
import csv
# ## Simonetti - Salles Da Cunha - Lucena Constraints model
# \begin{equation}
# \min \sum_{i \in V} x_i
# \end{equation}
# where $\mathscr{R}_0$ is the polyhedral region implied by:
#
# \begin{align}
#      \displaystyle \sum_{(i,j)\in E} y_{ij}  &=  \displaystyle \sum_{i \in V} x_i-1 \label{2a} \\
#     \displaystyle \sum_{(i,j) \in E(S)} y_{ij}  &\leq  \displaystyle
#           \sum_{i \in S\backslash \{j\}} x_i,   &S\subset V, j\in S \label{2b} \\
#      &\displaystyle \sum_{j\in \Gamma_i} x_j  \geq  1,  & i \in V \label{2c} \\
#     &y_{ij} \geq  0, & \forall(i,j)\in E \label{2d}\\
#      &0 \leq x_i  \leq  1, & \forall i\in V \label{2e}
# \end{align}


class Simonetti_SallesDaCunha_Lucena_Model:
    def __init__(self, V,E,A, status=True):
        self.V = V
        self.E = E
        self.A = A
        self.status = status
        self.model = Model("SSL")
        self.x = self.model.binary_var_dict(V, name='x')
        self.y = self.model.integer_var_dict(E, name='y')
        self.gamma_i = {}
        for i in V:
            self.gamma_i[i] = [i]
            for j in V:
                if i!=j and (i,j) in self.E or (j,i) in self.E:
                    self.gamma_i[i].append(j)
    def _build_model(self):
        self.model.minimize(self.model.sum(self.x[i] for i in self.V))
        # Constraints 1.3a
        self.model.add_constraint(self.model.sum(self.y)==self.model.sum(self.x)-1)
        # Constraints 1.3d
        self.model.add_constraints(self.y[i,j]>=0 for i,j in self.E)

        # Constraints 1.3e
        self.model.add_constraints(self.x[i]<=1 for i in self.V)
        self.model.add_constraints(self.x[i]>=0 for i in self.V)

        # Constraints 1.4a
        for v in self.V:
            self.model.add_constraint(self.model.sum(self.x[k] for k in self.gamma_i[v])-
                    self.model.sum(self.y[i,j] for i,j in self.E if i in self.gamma_i[v]
                    and j in self.gamma_i[v])>=1)

        # Constraints 1.3f
        self.model.add_constraints(self.y[i,j]<=self.x[i] for i,j in self.E)
        self.model.add_constraints(self.y[i,j]<=self.x[j] for i,j in self.E)

    def _update_constraints(self):
        #find if connected
        graph = {}
        self.active_vertices = [i for i in self.V if self.x[i].solution_value>0.9]
        self.active_edges = [(i,j) for i,j in self.E if self.y[i,j].solution_value>0.9]
        for i in self.active_vertices:
            graph[i]=[]
            for (j,k) in self.active_edges:
                if i==j:
                    graph[i].append(k)
                elif i==k:
                    graph[i].append(j)
        connected = np.zeros(len(self.active_vertices))
        i = self.active_vertices[0]
        connected[0] = 1
        for j in range(1,len(self.active_vertices)):
            #If exist a path from i to j
            if exist_path(graph, i, self.active_vertices[j]):
                connected[j] = 1
        #if connected, return True, optimal solution found
        if connected.all():
            return True
        cycles = []
        for node in self.active_vertices:
            cyclenodes = [[node]+path for path in dfs(graph, node, node)]
            for cycle in cyclenodes:
                if len(cycle)>0:
                    if (cycle not in cycles and cycle.reverse() not in cycles):
                        cycles.append(cycle)
            del graph[node]
            for key in graph:
                if node in graph[key]:
                    graph[key].remove(node)
        #else: add constraints to the model and return false
        if len(cycles)>0:
            for cycle in cycles:
                self.model.add_constraints(self.model.sum(self.y[i,k] for (i,k) in self.E
                if i in cycle and k in cycle)<=self.model.sum(self.x[i] for i in cycle if i!=j)
                for j in cycle)

        return False

    def solve_model(self):
        print("SSL")
        found_optimal = False
        self.iteration = 0
        start = time()*1000
        while not found_optimal:
            self.iteration+=1
            print("Iteration",self.iteration)
            self.model.parameters.timelimit = 3600 #No more than an hour
            self.model.parameters.mip.tolerances.mipgap = 0.05
            res = self.model.solve(clean_before_solve=True, log_output=self.status)

            found_optimal=self._update_constraints()
            if self.iteration>200 or time()*1000-start>3600000:
                break
        end =  time()*1000
        elapsed = int(round(end-start))
        self.write_info(elapsed, res)
        #print(self.model.objective_value)
        return res, self.active_vertices, self.active_edges

    def write_info(self, time, res):
        density = int(len(self.E)*2/(len(self.V)*len(self.V)-1)*100)
        filename = "results/SSL_"+str(len(self.V))+"_"+str(density)+".csv"
        with open(filename, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['SSL', len(self.V), len(self.E), time, self.model.objective_value, self.model.number_of_variables, self.model.number_of_constraints ])
        csvfile.close()


def exist_path(graph, start, end):
    fringe = [(start, [])]
    while fringe:
        state, path = fringe.pop()
        if path and state == end:
            return True
        for nextStep in graph[state]:
            if nextStep in path:
                continue
            fringe.append((nextStep, path+[nextStep]))
    return  False

def dfs(graph, start, end):
    fringe = [(start, [])]
    while fringe:
        state, path = fringe.pop()
        if path and state == end:
            if len(path)>3:
                yield path
            continue
        for next_state in graph[state]:
            if next_state in path:
                continue
            fringe.append((next_state, path+[next_state]))


def Simonetti_SallesDaCunha_Lucena(V,E,A, status=True):
    instance = Simonetti_SallesDaCunha_Lucena_Model(V,E,A,status)
    instance._build_model()
    return instance.solve_model()


#end
