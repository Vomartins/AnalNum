import numpy as np
import sys

class Cholesky:
    
    def __init__(self, A):
        self.__A = A
        
    def cholesky(self):
        #Checa se a matriz é quadrada.
        if(np.shape(self.__A)[0] != np.shape(self.__A)[1]):
            sys.stderr.write("ERRO: a matriz não é quadrada.")
            input()
            sys.exit()

        R = np.zeros((np.shape(self.__A)[0], np.shape(self.__A)[1])) #Fator de Cholesky

        for i in range(np.shape(self.__A)[0]):
            a = self.__A[i,i]
            if(i>0):
                a = a - np.dot(R[0:i,i], R[0:i,i])
            if(a<=0): #Verifica se é positiva definida.
                sys.stderr.write("ERRO: a matriz não é positiva definida.")
                input()
                sys.exit()
            else:
                R[i,i] = a**0.5

            if(i < np.shape(self.__A)[0]-1):
                for j in range(i+1, np.shape(self.__A)[1]):
                    R[i,j] = (self.__A[i,j] - np.dot(R[0:i,i],R[0:i,j]))/R[i,i]

        return R
    