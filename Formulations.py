from docplex.mp.model import Model


# ## Minimum dominating set function
#
#
# \min \sum_{i}^V x_i
#
# \sum_j^V a_{ij} x_j >= 1,  \forall i \in V

def MDS(V,E,A, status=True):
    mdl = Model('MDS')

    x = mdl.binary_var_dict(V, name='x')
    mdl.minimize(mdl.sum(x[i] for i in V))
    mdl.add_constraints(mdl.sum(A[i-1][j-1]*x[j] for j in V)>=1 for i in V)
    try:
        solution = mdl.solve(log_output=status)
        print(solution.solve_status)
        print(solution)
        active_vertices = [i for i in V if x[i].solution_value>0.9]
        density = int(len(E)*2/(len(V)*len(V)-1)*100)
        filename = "../results/MDS_"+str(len(V))+"_"+str(density)+".txt"
        with open(filename, 'a') as f:
            f.write(str(len(V))+" "+str(len(E))+"\n")
            f.write("\nObjective value: "+str(mdl.objective_value)+"\n")
            f.write(str(mdl.get_statistics()))
            f.write("\n---------------------------------\n\n")
        f.close()
    except AttributeError:
        # In case no solution exists
        active_vertices = []
    return solution, active_vertices
