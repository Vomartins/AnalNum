import Matriz
import numpy as np
from math import factorial as fact #Calculo de invH

#Cálculo do número de condicionamento analítico pela norma2(k2)
def cond(A, invA):
    return Matriz(A).norma2()*Matriz(invA).norma2

#Estimative inferior para o número de condicionamento pela norma2(k2)
def condest(A):
    n = np.shape(A)[0]
    a = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            a[i,j] = Vetor(A[:,i]).norma2()/Vetor(A[:,j]).norma2()
    return np.max(a)


#Inversa "analítica"
invH = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        invH[i,j] = ((-1)**(i+j))*fact(n+i)*fact(n+j)/(fact(n-j-1)*fact(n-i-1)*(i+j+1)*((fact(i)*fact(j))**2))