from cholesky import Cholesky
import numpy as np

def HilbertMatrix(n):
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = 1/(i+j+1)
    return A

A = HilbertMatrix(3)
R = Cholesky(A)

print('\nA =\n', A)
print('\nR = \n', R.cholesky())
print('\nR^t.R = \n', np.dot(R.cholesky().transpose(),R.cholesky()))

input()