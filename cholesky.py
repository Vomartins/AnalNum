import numpy as np
import sys

def Cholesky(A):
    #Checa se a matriz é quadrada.
    if(np.shape(A)[0] != np.shape(A)[1]):
        sys.stderr.write("ERRO: a matriz não é quadrada.")
        input()
        sys.exit()

    R = np.zeros((np.shape(A)[0], np.shape(A)[1])) #Fator de Cholesky
        
    for i in range(np.shape(A)[0]):
        a = A[i,i]
        if(i>0):
            a = a - np.dot(R[0:i,i], R[0:i,i])
        if(a<=0): #Verifica se é positiva definida.
            sys.stderr.write("ERRO: a matriz não é positiva definida.")
            input()
            sys.exit()
        else:
            R[i,i] = a**0.5

        if(i < np.shape(A)[0]-1):
            for j in range(i+1, np.shape(A)[1]):
                R[i,j] = (A[i,j] - np.dot(R[0:i,i],R[0:i,j]))/R[i,i]

    return R
     

def HilbertMatrix(n):
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = 1/(i+j+1)
    return A

def main():
    A = HilbertMatrix(3)
    R = Cholesky(A)

    print('\nA =\n', A)
    print('\nR = \n', R)
    print('\nR^t.R = \n', np.dot(R.transpose(),R))
   
    input()

main()