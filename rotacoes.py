import numpy as np
import sys

def QRrot(A):
    m = np.shape(A)[0]
    n = np.shape(A)[1]

    if(m<n):
        sys.stderr.write("ERRO: a matriz não admite fatoração QR.")
        input()
        sys.exit()

    Q = np.identity(m, dtype=float) #Q^t
    R = A.astype('float')

    for i in range(n):
        for j in range(i+1, m):
            if(R[j,i] == 0): #A matriz de rotação é a identidade.
                continue
            
            #Cálculo do cos() e sin()
            N = max([abs(R[i,i]),abs(R[j,i])]) #Fator de Normalização
            a = R[i,i]/N
            b = R[j,i]/N
            v = (a**2 + b**2)**0.5

            c = a/v
            s = b/v

            #Produto das matrizes
            #Matriz R
            Li = c*R[i] + s*R[j]
            Lj = c*R[j] - s*R[i]
            R[i] = Li.copy()
            R[j] = Lj.copy()

            #Matriz Q
            Li = c*Q[i] + s*Q[j]
            Lj = c*Q[j] - s*Q[i]
            Q[i] = Li.copy()
            Q[j] = Lj.copy()
            
    return Q.transpose(), R

def HilbertMatrix(n):
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = 1/(i+j+1)
    return A

def main():
    A = HilbertMatrix(13)
    Q,R = QRrot(A)

    B = np.dot(Q,R)
    B = A-B
    print(np.linalg.eigvals(np.dot(B.transpose(),B)).max())
   
    input()

main()

