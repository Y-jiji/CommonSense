import os

import gurobipy
import cvxpy as cp

env_g = gurobipy.Env()
env_g.setParam("OutputFlag", 0)
env_g.setParam("Threads", 1)

""" Solve x with min ||x||_1 s.t. Ax = y """
def compressive_solve(A, y):
    x = cp.Variable(A.shape[1])
    prob = cp.Problem(cp.Minimize(cp.norm(x, 1)),
                      [A @ x == y])
    prob.solve(solver=cp.GUROBI, env=env_g, verbose=False)
    return x.value

if __name__ == "__main__":
    import numpy as np
    A = np.array([[1, 2, 7], [3, 4, 1]])
    y = np.array([5, 6])
    print(compressive_solve(A, y))