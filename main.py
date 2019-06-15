import numpy as np
from LinSolve import LinSolve
from Vetor import Vetor, Dot
from Matriz import Matriz
from Fatoracao import Fat

n = 3

#Matriz de Hilbert
H = np.fromfunction(lambda i, j: 1/(i+j+1), (n,n), dtype=int)

#Inicialização do objeto Fat(), a partir dele podemos calcular as três fatorações da matriz de Hilbert.
A = Fat(H)

#Cálculo das fatorações a partir do objeto A e inicialização dos objetos Matriz com as fatorações obtidas.
R = Matriz(A.cholesky()) #Cholesky
q1, r1 = A.GS() #QR com Gram-Schmidt
Q1 = Matriz(q1)
R1 = Matriz(r1)
q2, r2 = A.rot() #QR rotação
Q2 = Matriz(q2)
R2 = Matriz(r2)

#Resolução do sistema HX=I com cholesky.
y1 = LinSolve(R.matriz.transpose(), I[:, 0]).triinf()
y2 = LinSolve(R.matriz.transpose(), I[:, 1]).triinf()
y3 = LinSolve(R.matriz.transpose(), I[:, 2]).triinf()

Y = np.array([y1, y2, y3])

x1 = LinSolve(R.matriz, Y[0, :]).trisup()
x2 = LinSolve(R.matriz, Y[1, :]).trisup()
x3 = LinSolve(R.matriz, Y[2, :]).trisup()

X = np.array([x1, x2, x3]).transpose()
If = np.dot(H, X)

print(X)
print(If)

#Resolução do sistema HX=I com Gram-Schmidt.
x1 = LinSolve(R1.matriz, Q1.transpose[:, 0]).trisup()
x2 = LinSolve(R1.matriz, Q1.transpose[:, 1]).trisup()
x3 = LinSolve(R1.matriz, Q1.transpose[:, 2]).trisup()

X = np.array([x1, x2, x3]).transpose()
If = np.dot(H, X)

print(X)
print(If)

#Resolução do sistema HX=I com rotações.
x1 = LinSolve(R2.matriz, Q2.transpose[:, 0]).trisup()
x2 = LinSolve(R2.matriz, Q2.transpose[:, 1]).trisup()
x3 = LinSolve(R2.matriz, Q2.transpose[:, 2]).trisup()

X = np.array([x1, x2, x3]).transpose()
If = np.dot(H, X)

print(X)
print(If)


input()