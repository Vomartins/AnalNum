import numpy as np
from LinSolve import LinSolve
from Vetor import Vetor, Dot
from Matriz import Matriz
from cholesky import Cholesky
from QR import QR

n = 3

A = np.fromfunction(lambda i, j: 1/(i+j+1), (n,n), dtype=int)
R = Cholesky(A)

print('\nA =\n', A)
print('\nR = \n', R.cholesky())
print('\nR^t.R = \n', np.dot(R.cholesky().transpose(),R.cholesky()))

input()