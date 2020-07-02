from docplex.mp.model import Model
from time import time
import csv

# ## Martin Constraints
# \begin{align}
#     \min &\displaystyle \sum_{i\in V} x_i\\
#     s.t. & \displaystyle \sum_{j\in V} a_{ij}x_j \geq 1 \\
#     x_i &\in \{0,1\}, & \forall i \in V
# \end{align}
# The Martin connectivity constraints are:
# \begin{align}
#     \displaystyle \sum_{(i,j)\in E}y_{ij} =\displaystyle \sum_{i\in V} x_i -1 \label{3a}\\
#     y_{ij}\leq x_i, y_{ij}\leq x_j, && \forall(i,j)\in E \label{3b}\\
#     z_{ij}^{k}\leq y_{ij},  z_{ij}^{k} \leq x_k, && \forall (i,j)\in E, k\in V \label{3c}\\
#     z_{ji}^{k}\leq y_{ij},   z_{ji}^{k} \leq x_k, && \forall(i,j)\in E, k \in V \label{3d}\\
#     y_{ij}- M(3-x_i-x_j-x_k) \leq& z_{ij}^{k}+ z_{ji}^{k} \leq y_{ij} +
#                               M(3-x_i-x_j-x_k), &\forall i,j,k \in V \label{3e}\\
#     1- M(2-x_i-x_j)  \leq &\displaystyle \sum_{k'\in V\backslash\{i,j\}}z_{ik'}^{j}+y_{ij}
#                   \leq 1+M(2-x_i-x_j), &\forall i,j \in V \label{3f}\\
#     y_{ij}, z_{ij}^{k}\in \{0,1\}, & \forall(i,j)\in E, k\in V, y_{ij} = 0, z_{ij}^{k}=0,
#                   \forall i,j,k \in V, &(i,j)\notin E \label{3g}
# \end{align}

# In[73]:
M = 10
class Martin_Model:
    def __init__(self, V,E,A, status=True):
        self.V = V
        self.E = E
        self.edges = [(i,j) for i in V for j in V]
        self.EE = []
        for i,j in E:
            self.EE.append((i,j))
            self.EE.append((j,i))
        self.A = A
        self.Z = [(i,j,k) for i in V for j in V for k in V]
        self.status = status
        self.model = Model("Martin")
        self.x = self.model.binary_var_dict(V, name="x")
        self.y = self.model.binary_var_dict(self.edges, name="y")
        self.z = self.model.binary_var_dict(self.Z, name="z")

    def _build_model(self):
        #Objective function
        self.model.minimize(self.model.sum(self.x[i] for i in self.V))
        self.model.add_constraints(self.model.sum(self.A[i-1][j-1]*self.x[j]
                                for j in self.V)>=1 for i in self.V)
        #Constraint 3a
        self.model.add_constraint(self.model.sum(self.y[i,j] for i,j in self.E)
                                == self.model.sum(self.x[i] for i in self.V)-1)
        #Constraint 3b
        self.model.add_constraints(self.y[i,j]<=self.x[i] for i,j in self.E)
        self.model.add_constraints(self.y[i,j]<=self.x[j] for i,j in self.E)
        #Constraint 3c
        self.model.add_constraints(self.z[i,j,k]<=self.y[i,j] for i,j in self.E for k in self.V)
        self.model.add_constraints(self.z[i,j,k]<=self.x[k] for i,j in self.E for k in self.V)
        #Constraint 3d
        self.model.add_constraints(self.z[j,i,k]<=self.y[i,j] for i,j in self.E for k in self.V)
        self.model.add_constraints(self.z[j,i,k]<=self.x[k] for i,j in self.E for k in self.V)
        #Constraint 3e
        self.model.add_constraints(self.y[i,j]-M*(3-self.x[i]-self.x[j]-self.x[k])
                    <=self.z[i,j,k]+self.z[j,i,k] for i in self.V for j in self.V for k in self.V)
        self.model.add_constraints(self.z[i,j,k]+self.z[j,i,k] <= self.y[i,j]+
                    M*(3-self.x[i]-self.x[j]-self.x[k])  for i in self.V
                    for j in self.V for k in self.V )
        #Constraint 3f
        self.model.add_constraints(1-M*(2-self.x[i]-self.x[j])<=self.model.sum(self.z[i,k,j]
                    for k in self.V if k!=i and k!=j)+self.y[i,j] for i in self.V for j in self.V )
        self.model.add_constraints(1+M*(2-self.x[i]-self.x[j])>=self.model.sum(self.z[i,k,j]
                    for k in self.V if k!=i and k!=j)+self.y[i,j] for i in self.V for j in self.V )
        #Constraint 3g
        self.model.add_constraints(self.y[i,j]==0 for i in self.V for j in self.V
                                if (i,j) not in self.E and (j,i) not in self.E)
        self.model.add_constraints(self.z[i,j,k]==0 for i in self.V for j in self.V
                    for k in self.V if  (i,j) not in self.EE)

    def solve_model(self):
        print("Martin")
        self.model.parameters.timelimit = 3600
        #self.model.parameters.mip.tolerances.mipgap = 0.5
        start = time()*1000
        res = self.model.solve(clean_before_solve=True, log_output=self.status)
        end = time()*1000
        elapsed = int(round(end-start))
        self.write_info(elapsed, res)
        active_vertices = [i for i in self.V if self.x[i].solution_value>0.9]
        active_edges = [(i,j) for i,j in self.E if self.y[i,j].solution_value>0.9]
        return res, active_vertices, active_edges
    def write_info(self, time, res):
        density = int(len(self.E)*2/(len(self.V)*len(self.V)-1)*100)
        filename = "../results/Martin_"+str(len(self.V))+"_"+str(density)+".csv"
        with open(filename, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['SSL_L', len(self.V), len(self.E), time, self.model.objective_value, self.model.number_of_variables, self.model.number_of_constraints ])
        csvfile.close()


def Martin(V,E,A, status=True):
    instance = Martin_Model(V,E,A,status)
    instance._build_model()
    try:
        return instance.solve_model()
    except AttributeError or DOcplexException:
        print("No solution found")
        active_vertices = []
        active_edges = []
    return solution, active_vertices, active_edges
