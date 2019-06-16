import numpy as np
import sys
from Vetor import Dot, Vetor

class Fat:

    def __init__(self,A):
        self.__A = A
        
    #Fatoração de Cholesky.   
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
    
    #Fatoração QR com Gram-Schmidt.
    def GS(self):
        m, n = self.__A.shape
        Q = np.zeros((m, n))
        R = np.zeros((n, n))

        for k in range(n):
            v = self.__A[:, k]
            
            v_aux = np.zeros(m) #Soma das projeções de v sobre os k-1 vetores ortonormais

            for i in range(k):
                q = Q[:, i]
                R[i, k] = Dot(q, v).dot()
                v_aux = v_aux + R[i, k] * q
        
            v = v - v_aux

            norma = Vetor(v).norma2()
           
            Q[:, k] = v/norma
            R[k, k] = norma
        return Q, R


    #Fatoração QR com rotação.
    def rot(self):
        m = np.shape(self.__A)[0] #Número de linhas de A
        n = np.shape(self.__A)[1] #Número de colunas de A

        if(m<n):
            sys.stderr.write("ERRO: a matriz não admite fatoração QR.")
            input()
            sys.exit()

        Q = np.identity(m, dtype=float) #Q^t
        R = self.__A.astype('float')

        for i in range(n):
            for j in range(i+1, m):
                if(R[j,i] == 0): #A matriz de rotação é a identidade.
                    continue
                
                #Cálculo do cos() e sin()
                N = max([abs(R[i,i]),abs(R[j,i])]) #Fator de Normalização
                a = R[i,i]/N
                b = R[j,i]/N
                v = (a**2 + b**2)**0.5

                c = a/v #cos
                s = b/v #sin

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