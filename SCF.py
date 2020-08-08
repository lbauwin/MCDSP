from docplex.mp.model import Model
from time import time
import csv

# ## Single Commodity Flow Constraints
# \begin{align}
#     \displaystyle \sum_{i\in V} r_i &=1 \label{4a}\\
#     r_i &\leq x_i, &\forall i \in V \label{4b}\\
#     f_{ij}&\geq 0, & \forall(i,j) \in E \cup E' \label{4c}\\
#     f_{ij}\leq x_i \displaystyle \sum_{k\in V} x_k&,  f_{ij}\leq
#           \displaystyle x_j\sum_{k\in V} x_k, &\forall(i,j) \in E\cup E' \label{4d}\\
#     \displaystyle \sum_{j} f_{ji}&\leq n(1-r_i), &\forall i\in V \label{4e}\\
#     \displaystyle \sum_{j}f_{ji} - \displaystyle \sum_{j}f_{ij} &=x_i-r_i
#           \displaystyle \sum_{j\in V} x_j, & \forall i\in V \label{4f}\\
#     r_i &\in \{0,1\}, & \forall i\in V \label{4g}
# \end{align}

class Single_Commodity_Flow_Model:
    def __init__(self, V,E,A, status=True):
        self.V = V
        self.E = E
        self.A = A
        self.edges=[]
        self.status = status
        for i in range(len(A)):
            for j in range(len(A[i])):
                if A[i][j]==1 and i!=j:
                    self.edges.append((i+1,j+1))
        self.model = Model("SCF")
        self.x = self.model.binary_var_dict(V, name="x")
        self.r = self.model.binary_var_dict(V, name='r')
        self.f = self.model.integer_var_dict(self.edges, name='f')

    def _build_model(self):
        #Objective function
        self.model.minimize(self.model.sum(self.x[i] for i in self.V))
        self.model.add_constraints(self.model.sum(self.A[i-1][j-1]*self.x[j]
                                for j in self.V)>=1 for i in self.V)
        #Constraint 1.10a
        self.model.add_constraint(self.model.sum(self.r)==1)
        #Constraint 1.10b
        self.model.add_constraints(self.r[i]<=self.x[i] for i in self.V)
        #Constraint 1.10c
        self.model.add_constraints(self.f[i,j]>=0 for i,j in self.edges)
        #Constraint 1.10d
        self.model.add_indicator_constraints(self.model.indicator_constraint(self.x[i],self.f[i,j]
                    <=self.model.sum(self.x[k] for k in self.V),1) for i,j in self.edges)
        self.model.add_indicator_constraints(self.model.indicator_constraint(self.x[i],self.f[i,j]
                    ==0,0)for i,j in self.edges)
        self.model.add_indicator_constraints(self.model.indicator_constraint(self.x[j],self.f[i,j]
                    <=self.model.sum(self.x[k] for k in self.V),1) for i,j in self.edges)
        self.model.add_indicator_constraints(self.model.indicator_constraint(self.x[j],self.f[i,j]
                    ==0, 0) for i,j in self.edges)
        #Constraint 1.10e
        self.model.add_constraints(self.model.sum(self.f[j,i] for j in self.V if (j,i) in self.edges)
                    <=len(self.V)*(1-self.r[i]) for i in self.V)
        #Constraint 1.10f
        self.model.add_indicator_constraints(self.model.indicator_constraint(self.r[i],
                    self.model.sum(self.f[j,i] for j in self.V if (j,i) in self.edges)-
                    self.model.sum(self.f[i,j] for j in self.V if (i,j) in self.edges)
                    ==self.x[i]-self.model.sum(self.x[j] for j in self.V),1)for i in self.V)
        self.model.add_indicator_constraints(self.model.indicator_constraint(self.r[i],
                    self.model.sum(self.f[j,i] for j in self.V if (j,i) in self.edges)-
                    self.model.sum(self.f[i,j] for j in self.V if (i,j) in self.edges)
                    ==self.x[i],0) for i in self.V)

    def solve_model(self):
        print("SCF")
        self.model.parameters.timelimit = 3600
        self.model.parameters.mip.tolerances.mipgap = 0.05
        start = time()*1000
        res = self.model.solve(clean_before_solve=True, log_output=self.status)
        end = time()*1000
        #print(self.model.objective_value)

        elapsed = int(round(end-start))
        self.write_info(elapsed, res)
        active_vertices = [i for i in self.V if self.x[i].solution_value>0.9]
        active_edges = [(i,j) for i,j in self.E if self.f[i,j].solution_value>0.9 or self.f[j,i].solution_value>0.9]
        return res, active_vertices, active_edges
        
    def write_info(self, time, res):
        density = int(len(self.E)*2/(len(self.V)*len(self.V)-1)*100)
        filename = "../results/SCF_"+str(len(self.V))+"_"+str(density)+".csv"
        with open(filename, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['SCF', len(self.V), len(self.E), time, self.model.objective_value, self.model.number_of_variables, self.model.number_of_constraints ])
        csvfile.close()

def Single_Commodity_Flow(V,E,A, status=True):
    instance = Single_Commodity_Flow_Model(V,E,A,status)
    instance._build_model()
    try:
        return instance.solve_model()
    except AttributeError or DOcplexException:
        print("No solution found")
        active_vertices = []
        active_edges = []
    return solution, active_vertices, active_edges
