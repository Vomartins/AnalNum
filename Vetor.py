import numpy as np

class Dot:
    
    def __init__(self, x, y):
        self.__x = Vetor(x)
        self.__y = Vetor(y)
    
    def dot(self):
        n = np.size(self.__x.vetor)
        soma = 0
        for i in range(n):
            soma = soma + self.__x.vetor[i]*self.__y.vetor[i]
        return soma

class Vetor:
    
    def __init__(self, X):
        self.__X = X
        
    @property
    def vetor(self):
        return self.__X
    
    @property
    def dim(self):
        return np.shape(self.__X)[0]
        
    #Norma 1 de vetor.
    def norma1(self):
        n = 0 
        mod_x = np.abs(self.__X)
        for i in range(np.shape(self.__X)[0]):
            n = n + mod_x[i]
        return n
    
    #Norma 2 de vetor.
    def norma2(self):
        n = 0 
        b = abs(self.__X).max()
        if b == 0:
            n = 0
        else:
            x = (1/b)*self.__X
            for i in range(np.shape(self.__X)[0]):
                n = n + x[i]**2
            n = b*n**(1/2)
        return n