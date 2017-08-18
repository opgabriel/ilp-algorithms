# Algoritmo primal-dual para o problema de set-cover por Gabriel Oliveira

import numpy as np
import sys

# Recebe vetor c e matriz A da entrada
A = np.loadtxt(sys.argv[1], dtype = float, skiprows = 2)

# Extrai vetor c e matriz A como submatrizes
cT = np.asmatrix(A[0][:])
A = np.asmatrix(A[1:][:])

# Extrai dimensões da matriz
(nLinhas, nColunas) = A.shape

# Inicializa vetores x, y e vetor de 1s
x = np.asmatrix(np.zeros(nLinhas)).getT()     # Solução sempre viavel da dual
y = np.asmatrix(np.zeros(nColunas)).getT()    # Só será viável para a primal na última iteração
vetUm = np.asmatrix(np.ones(nLinhas)).getT() 

# DUAL: Inicializa A^T e c
c = cT.getT()
AT = A.getT()

# Loop do algoritmo
Ay = np.dot(A, y)
while not np.greater_equal(Ay, vetUm).all():	# Enquanto houver elemento descoberto (A*y < 1)
	naoCoberto = Ay.argmin()   # Recebe indice do primeiro elemento que nao esta coberto -- Ay[naoCoberto] = 0

	# Aumenta x_a ao maximo de forma que AT*x <= c
	ATx = np.dot(AT, x)
	while np.less_equal(ATx, c).all():
		x[naoCoberto] += 0.001 # Aumenta ate achar igualdade
		ATx = np.dot(AT, x)

		# Verifica se é a última iteração
		x_temp = x.copy()
		x_temp[naoCoberto] += 1
		ATx_temp = np.dot(AT, x_temp)
		if np.greater(ATx_temp, c).any(): # Se alguma linha de A^T*x > c
			linhaIgualdade = np.greater(ATx_temp, c).tolist().index([True]) # Linha onde ocorreu a igualdade
			break
	
	
	y[linhaIgualdade] = 1 	# Adiciona o subconjunto à resposta

	Ay = np.dot(A, y)

	# Impressão dos vetores intermediários
	print("y = ", np.around(np.array(y.T)[0].tolist(), decimals = 1))
	print("x = ", np.around(np.array(x.T)[0].tolist(), decimals = 1), "\n\n")

# Calculo do custo da cobertura
custo = np.dot(cT, y)
print("Custo da cobertura: " + str(np.around(float(custo), decimals = 1)))
