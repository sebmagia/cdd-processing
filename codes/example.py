import numpy as np

from scipy.sparse.linalg import LinearOperator,aslinearoperator

def mv(v):
    return np.array([2*v[0],3*v[1]])

A=LinearOperator((2,2),matvec=mv)
print(A)
