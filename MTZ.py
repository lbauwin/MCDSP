from docplex.mp.model import Model
from time import time
import csv
# ### Miller Tucker Zemlin Constraints
#
# \begin{align}
#     \displaystyle \sum_{i\in V} y_{n+2,i} &=1 \label{5a}\\
#     \displaystyle \sum_{i:(i,j)\in A} y_{ij} &=1, & \forall j\in V \label{5b}\\
#     y_{n+1,i}+y_{i,j} & \leq 1, & \forall (i,j)\in E\cup E' \label{5c}\\
#     (n+1)y_{ij} + u_i-u_j+(n-1)y_{ji} &\leq n, & \forall (i,j)\in E\cup E' \label{5d}\\
#     (n+1)y_{ij} + u_i-u_j & \leq n, & \forall(i,j)\in A\setminus (E\cup E') \label{5e}\\
#     y_{n+1,n+2} &=1 \label{5f}\\
#     u_{n+1} &=0 \label{5g}\\
#     1 \leq u_i &\leq n+1, & i\in V \cup \{n+2\} \label{5h}\\
#     x_i &= 1-y_{n+1,i}, & \forall i\in V \label{5i}
# \end{align}
class Miller_Tucker_Zemlin_Model:
    def __init__(self, V,E,A, status=True):
        self.V = V
        self.E = E
        self.A = A
        self.status = status
        self.model = Model("MTZ")
        self.v=len(V)
        self.U = V+[self.v+1]+[self.v+2]
        self.VV = V+[self.v+2]
        self.YY = []
        for i in V:
            self.YY.append((self.v+1,i))
            self.YY.append((self.v+2,i))
        self.YY = self.YY+[(self.v+1,self.v+2)]
        self.EE = []
        for i,j in E:
            self.EE.append((i,j))
            self.EE.append((j,i))
        self.Y = self.EE+self.YY
        # Variables
        self.x = self.model.binary_var_dict(V, name='x')
        self.y = self.model.binary_var_dict(self.Y, name='y')
        self.u = self.model.integer_var_dict(self.U, name='u')

    def _build_model(self):
        self.model.minimize(self.model.sum(self.x[i] for i in self.V))
        self.model.add_constraints(self.model.sum(self.A[i-1][j-1]*self.x[j]
                    for j in self.V)>=1 for i in self.V)
        #Constraint 1.12a
        self.model.add_constraint(self.model.sum(self.y[self.v+2,i] for i in self.V)==1)
        #Constraint 1.12b
        self.model.add_constraints(self.model.sum(self.y[i,j] for i in self.U if (i,j) in self.Y)
                    ==1 for j in self.V)
        #Constraint 1.12c
        self.model.add_constraints(self.y[self.v+1,i]+self.y[i,j]<=1 for i,j in self.EE)
        #Constraint 1.12d
        self.model.add_constraints((self.v+1)*self.y[i,j]+self.u[i]-self.u[j]+(self.v-1)*self.y[j,i]
                    <=self.v for i,j in self.EE)
        #Constraint 1.12e
        self.model.add_constraints((self.v+1)*self.y[i,j]+self.u[i]-self.u[j]<=self.v for i,j in self.YY)
        #Constraint 1.12f
        self.model.add_constraint(self.y[self.v+1,self.v+2]==1)
        #Constraint 1.12g
        self.model.add_constraint(self.u[self.v+1]==0)
        #Constraint 1.12i
        self.model.add_constraints(1<=self.u[i] for i in self.VV)
        self.model.add_constraints(self.u[i]<=self.v+1 for i in self.VV)
        #Constraint 1.12i
        self.model.add_constraints(self.x[i]==1-self.y[self.v+1,i] for i in self.V)


    def solve_model(self):
        print("MTZ")
        self.model.parameters.timelimit = 3600
        self.model.parameters.mip.tolerances.mipgap = 0.05
        start = time()*1000
        res = self.model.solve(clean_before_solve=True, log_output=self.status)
        end = time()*1000
        if res == None:
            print('infeasible')
            return
        print(self.model.objective_value)

        elapsed = int(round(end-start))
        self.write_info(elapsed, res)
        active_vertices = [i for i in self.V if self.x[i].solution_value>0.9]
        active_edges = [(i,j) for i,j in self.E if self.y[i,j].solution_value>0.9]
        return res, active_vertices, active_edges

    def write_info(self, time, res):
        density = int(len(self.E)*2/(len(self.V)*len(self.V)-1)*100)
        filename = "results/MTZ_"+str(len(self.V))+"_"+str(density)+".csv"
        with open(filename, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['MTZ', len(self.V), len(self.E), time, self.model.objective_value, self.model.number_of_variables, self.model.number_of_constraints ])
        csvfile.close()


def Miller_Tucker_Zemlin(V,E,A, status=True):
    instance = Miller_Tucker_Zemlin_Model(V,E,A,status)
    instance._build_model()
    try:
        return instance.solve_model()
    except AttributeError or DOcplexException:
        print("No solution found")
        active_vertices = []
        active_edges = []
    return solution, active_vertices, active_edges
