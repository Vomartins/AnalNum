import numpy as np

class Matriz:
    
    def __init__(self, X):
        self.__X = X
        
    @property
    def matriz(self):
        return self.__X
    
    @property
    def dim(self):
        return np.shape(self.__X)
    
    @property
    def transpose(self):
        return self.__X.transpose()
    
    #Norma 1 de matriz.
    def norma1(self):
        N = np.zeros((np.shape(self.__X)[0], 1))
        for i in range(np.shape(self.__X)[0]):
            col = np.abs(self.__X[:,i])
            for j in range(np.shape(self.__X)[0]):
                N[i, 0] = N[i, 0] + col[j]
        n = N.max()
        return n
    
    #Norma 2 de matriz.
    def norma2(self):
        return (np.linalg.eigvals(self.__X.transpose()*self.__X).max())**(1/2)
    