import os

print('\n')
print('*************** Análise  Numérica ***************')
print('Você pode:')
print('(1) Calcular a inversa de uma matriz de Hilbert de dimensão n \n ou \n(2) Calcular as constantes de condicionamento de uma matriz de Hilbert de ordem n')
res = int(input('O que deseja fazer?'))

if res == 1:
    os.system("inversa.py")
elif res == 2:
    os.system("condicionamento.py")
else:
    print('Hmm, você não tem essa opção :( .')