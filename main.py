import numpy as np
from LinSolve import LinSolve
from Matriz import Matriz
from Vetor import Dot, Vetor
from Fatoracao import Fat
from math import factorial as fact #Calculo de invH

print('\n')
print('*************** Análise  Numérica ***************')
n = int(input('Qual a dimensão da matriz de Hilbert desejada? '))
print('Códigos dos Algortimos:')
print('    (1) Cholesky \n    (2) QR com Gram-Schmidt \n    (3) QR com rotação')
alg = int(input('Qual o algoritmo desejado? '))
print('\n')

w = []
for i in range(n):
    if (i%2 == 0):
        w.append(0)
    else:
        w.append(1)
w = np.array(w)

#Cálculo do número de condicionamento analítico pela norma1(k1)
def cond(A, invA):
    print(Matriz(A).norma1())
    print(Matriz(invA).norma1())
    return Matriz(A).norma1()*Matriz(invA).norma1()

#Estimative inferior para o número de condicionamento pela norma1(k1)
def condest(A):
    n = np.shape(A)[0]
    a = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            a[i,j] = Vetor(A[:,i]).norma1()/Vetor(A[:,j]).norma1()
    return np.max(a)

def condest2(A, invA, w):
    return Matriz(A).norma1()*Vetor(np.dot(invA,w)).norma1()/Vetor(w).norma1()

#Matriz de Hilbert
H = np.fromfunction(lambda i, j: 1/(i+j+1), (n,n), dtype=int)
#Inversa "analítica"
invH = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        invH[i,j] = ((-1)**(i+j))*fact(n+i)*fact(n+j)/(fact(n-j-1)*fact(n-i-1)*(i+j+1)*((fact(i)*fact(j))**2))

I = np.identity(n)        
                
#Inicialização do objeto Fat(), a partir dele podemos calcular as três fatorações da matriz de Hilbert.
A = Fat(H)
      
#Cálculo das fatorações de acordo com o algoritmo escolhido a partir do objeto A e inicialização dos objetos Matriz com as fatorações obtidas.
if alg == 1:
    R = Matriz(A.cholesky()) #Cholesky
    #Resolução do sistema HX=I com cholesky.
    Y = np.zeros((1, n))
    for i in range(n):
        y = LinSolve(R.matriz.transpose(), I[:, i]).Solve()
        Y = np.vstack((Y, y))
    Y = np.delete(Y, 0, 0)

    X = np.zeros((1, n))
    for i in range(n):
        x = LinSolve(R.matriz, Y[i, :]).Solve()
        X = np.vstack((X, x))
    X = np.delete(X, 0, 0)
    In = Matriz(np.dot(H, X) - I)
    r = In.norma1()
    r1 = Matriz(X-invH).norma1()
    
    k0 = np.around(cond(H, invH),decimals=3)
    k1 = np.around(cond(H, X), decimals=3)
    k2 = np.around(condest(H), decimals=3)
    k3 = np.around(condest2(H,X,w), decimals=3)
    
    if n <= 3:
        print('A fatoração de Cholesky é a seguinte:')
        print('Matriz R:')
        print(R.matriz)
        print('Matriz Rt:')
        print(R.transpose)
    if n < 6:
        print('Com a fatoração de Cholesky, a inversa da matriz de Hibert de ordem {} é:'.format(n))
        print(X)
        print('H*X - I:')
        print(In.matriz)
    
    print('Distância em relação à inversa exata: r1 = ', r1)    
    print('Raio da inversa: r = {}'.format(r))
    print('k0 = ',k0)
    print('k1 = ',k1)
    print('k2 = ',k2)
    print('k3 = ',k3)    
    
elif alg == 2:
    q1, r1 = A.GS() #QR com Gram-Schmidt
    Q1 = Matriz(q1)
    R1 = Matriz(r1)
    #Resolução do sistema HX=I com Gram-Schmidt.
    X = np.zeros((1,n))
    for i in range(n):
        x = LinSolve(R1.matriz, Q1.transpose[:, i]).Solve()
        X = np.vstack((X,x))
    X = np.delete(X, 0, 0)
    In = Matriz(np.dot(H, X) - I)
    r = In.norma1()
    r1 = Matriz(X-invH).norma1()
    
    k0 = np.around(cond(H, invH),decimals=3)
    k1 = np.around(cond(H, X), decimals=3)
    k2 = np.around(condest(H), decimals=3)
    k3 = np.around(condest2(H,X,w), decimals=3)
    
    if n <= 3:
        print('A fatoração de QR com Gram-Schmidt é a seguinte:')
        print('Matriz Q:')
        print(Q1.matriz)
        print('Matriz R:')
        print(R1.matriz)
    if n < 6:  
        print('Com a fatoração de QR com Gram-Schmidt, a inversa da matriz de Hibert de ordem {} é:'.format(n))
        print(X)
        print('H*X - I:')
        print(In.matriz)
    print('Raio da inversa: r = {}'.format(r))
    print('Distância em relação à inversa exata: r1 = ', r1)
    print('k0 = ',k0)
    print('k1 = ',k1)
    print('k2 = ',k2)
    print('k3 = ',k3)   
    
elif alg == 3:
    q2, r2 = A.rot() #QR rotação
    Q2 = Matriz(q2)
    R2 = Matriz(r2)
    #Resolução do sistema HX=I com rotações.
    X = np.zeros((1,n))
    for i in range(n):
        x = LinSolve(R2.matriz, Q2.transpose[:, i]).Solve()
        X = np.vstack((X,x))
    X = np.delete(X, 0, 0)
    In = Matriz(np.dot(H, X) - I)
    r = In.norma1()
    r1 = Matriz(X-invH).norma1()
    
    k0 = np.around(cond(H, invH),decimals=3)
    k1 = np.around(cond(H, X), decimals=3)
    k2 = np.around(condest(H), decimals=3)
    k3 = np.around(condest2(H,X,w), decimals=3)
    
    if n <= 3:
        print('A fatoração de QR com rotações é a seguinte:')
        print('Matriz Q:')
        print(Q2.matriz)
        print('Matriz R:')
        print(R2.matriz)
    if n < 6:
        print('Com a fatoração de QR com rotações, a inversa da matriz de Hibert de ordem {} é: '.format(n))
        print(X)
        print('H*X - I:')
        print(In.matriz)
    print('Raio da inversa: r = {}'.format(r))
    print('Distância em relação à inversa exata: r1 = ', r1)
    print('k0 = ',k0)
    print('k1 = ',k1)
    print('k2 = ',k2)
    print('k3 = ',k3)   
    
else:
      print('Hmm, não entendi :( .')
