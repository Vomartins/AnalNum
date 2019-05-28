import numpy as np

class LinSolve:
    
    def __init__(self, A, b):
        self.__A = A
        self.__b = b
        
    #Identifica qual o tipo da matriz e chama a função correta    
    def Solve(self):
        
        if self.__A[1,0] == 0 and self.__A[0,1] != 0:
            x1 = self.trisup()
            if (np.dot(self.__A, x1) - self.__b == np.zeros(np.shape(self.__b))).all():
                return x1
            else:
                return "A matriz não é triangular."
        elif self.__A[0,1] == 0 and self.__A[1,0] != 0:
            x2 = self.triinf()
            if (np.dot(self.__A, x2) - self.__b == np.zeros(np.shape(self.__b))).all():
                return x2
            else: 
                return "A matriz não é triangular."
        else:
            return "A matriz não é triangular."
        
    #Solução direta do sistem linear com matriz A triangular superior
    def trisup(self):
        x = np.zeros((np.shape(self.__A)[0],1))
        s = 0
        for j in range(np.shape(self.__A)[0]):
            x[np.shape(self.__A)[0]-(j+1), 0] = (self.__b[np.shape(self.__A)[0]-(j+1),0] - s)/ self.__A[np.shape(self.__A)[0]-(j+1),np.shape(self.__A)[0]-(j+1)]
            s = 0
            if j < np.shape(self.__A)[0]-1 :
                for i in range(j+1):
                    s = s + x[np.shape(self.__A)[0]-(i+1), 0]*self.__A[np.shape(self.__A)[0]-(i+2), np.shape(self.__A)[0]-(i+1)]
        return x

    #Solução direta do sistem linear com matriz A triangular inferior
    def triinf(self):
        x = np.zeros((np.shape(self.__A)[0],1))
        s = 0
        for j in range(np.shape(self.__A)[0]):
            x[j, 0] = (self.__b[j,0] - s)/ self.__A[j,j]
            s = 0
            if j < np.shape(self.__A)[0]-1:
                for i in range(j+1):
                    s = s + x[i, 0]*self.__A[j+1, i]
        return x