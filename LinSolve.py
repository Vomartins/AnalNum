import numpy as np

class LinSolve:
    
    def __init__(self, A, b):
        self.__A = A
        self.__b = b
        
    #Identifica qual o tipo da matriz e chama a função correta
    def Solve(self):
        S = 0
        for i in range(0, np.shape(self.__A)[0]):
            for j in range(i+1, np.shape(self.__A)[0]):
                S = S + self.__A[i, j]

        if S <= 10**(-10):
            return self.triinf()
        else:
            return self.trisup()
        '''
        if self.__A[1,0] == 0 and self.__A[0,1] != 0:
            x1 = self.trisup()
            if (np.dot(self.__A, x1) - self.__b <= np.full(np.shape(self.__b),10**(-10))).all():
                return x1
            else:
                return "A matriz não é triangular."
        elif self.__A[0,1] == 0 and self.__A[1,0] != 0:
            x2 = self.triinf()
            if (np.dot(self.__A, x2) - self.__b <= np.full(np.shape(self.__b),10**(-10))).all():
                return x2
            else: 
                return "A matriz não é triangular."
        else:
            return "A matriz não é triangular."
        '''
    #Solução direta do sistem linear com matriz A triangular superior
    def trisup(self):
        x = np.zeros((1,np.shape(self.__A)[0]))
        s = 0
        for j in range(np.shape(self.__A)[0]):
            x[0,np.shape(self.__A)[0]-(j+1)] = (self.__b[np.shape(self.__A)[0]-(j+1)] - s)/ self.__A[np.shape(self.__A)[0]-(j+1),np.shape(self.__A)[0]-(j+1)]
            s = 0
            if j < np.shape(self.__A)[0]-1:
                for i in range(j+1):
                    s = s + x[0,np.shape(self.__A)[0]-(i+1)]*self.__A[np.shape(self.__A)[0]-(j+2), np.shape(self.__A)[0]-(i+1)]
        return x

    #Solução direta do sistem linear com matriz A triangular inferior
    def triinf(self):
        x = np.zeros((1,np.shape(self.__A)[0]))
        s = 0
        for j in range(np.shape(self.__A)[0]):
            x[0,j] = (self.__b[j] - s)/ self.__A[j,j]
            s = 0
            if j < np.shape(self.__A)[0]-1:
                for i in range(j+1):
                    s = s + x[0,i]*self.__A[j+1, i]
        return x