""" Universidade Federal de Minas Gerais
	Instituto de Ciencias Exatas - ICEx
	Departamento de Ciencia da Computacao - DCC
	DCC035: Pesquisa Operacional - Prof. Gabriel Coutinho

   	Trabalho Pratico 2 - Programacoes Inteiras
   	ALGORITMO DE FORD-FULKERSON

   	Gabriel Pereira de Oliveira - 2015004160
"""

import numpy as np
import sys

# Função que atualiza o grafo residual
def atualizaGrafoResidual(G, f, c):
	nArestas = G.shape[1] 	# Numero de arestas é o numero de colunas da matriz de incidência
	nVertices = G.shape[0]	# Numero de vertices é o numero de linhas da matriz de incidência

	Gr = np.zeros((nVertices, nVertices))
	for aresta in range(nArestas):
		# Captura os indices dos vertices de chegada e saida da aresta
		u = int(np.where(G.getT()[aresta] == -1)[1])
		v = int(np.where(G.getT()[aresta] == 1)[1])

		# Cria matriz de adjacencia com os fluxos atualizados
		f_a = f.tolist()[0][aresta]
		c_a = c.tolist()[0][aresta]
		
		if f_a < c_a:
			Gr[u][v] = c_a - f_a
		if f_a > 0:
			Gr[v][u] = f_a

	return Gr


# Implementa de fato a visita DFS
def visita(G, visitados, s, t, maiorFluxoPossivel, caminho):
	achouCaminho = False

	# Neste caso, conseguiu chegar até t
	if s == t:
		achouCaminho = True


	# Marca o vértice atual como visitado
	visitados[s] = True
	caminho.append(s)
	#print("Visitando vertice ", s, "	|	Fluxo maximo: ", maiorFluxoPossivel % 1000, "	| 	Caminho ate agora: ", caminho)

	# Para cada vizinho de s, faz a visita (só procura se nao tiver achado caminho)
	if not achouCaminho:
		for x in range(len(visitados)):
			if (not visitados[x]) and (G[s][x] > 0):
				
				# Atualiza o valor do maior fluxo possivel
				if G[s][x] < maiorFluxoPossivel:
					maiorFluxoPossivel = G[s][x]

				(achouCaminho, fluxo, caminho) = visita(G, visitados, x, t, maiorFluxoPossivel, caminho)

				# Se tiver achado algum caminho, pode parar
				if achouCaminho:
					maiorFluxoPossivel = fluxo
					break

		if not achouCaminho:
			caminho.remove(s)

	return (achouCaminho, maiorFluxoPossivel, caminho)


# Função que acha um caminho entre s e t via busca por profundidade (DFS: depth-first search)
def achaCaminhoDFS(G, s, t):
	#print("\nIniciando busca de um st-caminho entre " + str(s) + " e " + str(t))
	visitados = np.full((1, G.shape[1]), False).tolist()[0]
	return visita(G, visitados, s, t, sys.maxsize, [])


# Função que, dado um caminho de vertices, encontra os indices das arestas correspondentes
def achaArestasCaminho(G, caminho):
	arestas = []
	G_aux = G.tolist()
	#print("achaArestasCaminho ---> Caminho de vertices recebido: ", caminho)
	
	# Pra cada vertice no caminho, pesquisa qual aresta o conecta com o prox vertice do caminho
	for vertice in caminho:
		ind_proximo = caminho.index(vertice) + 1

		# Caso o vertice seja o último do caminho
		if ind_proximo >= len(caminho):
			break

		proximo = caminho[ind_proximo]

		# Senão, procura qual aresta
		achouAresta = False
		for indiceAresta in range(G.shape[1]):
			if (G_aux[vertice][indiceAresta] == -1) and (G_aux[proximo][indiceAresta] == 1):
				#print("Adiciona aresta que liga o vertice ", vertice, " ao vertice ", proximo, ": Aresta ", indiceAresta)
				arestas.append(indiceAresta)
				achouAresta = True
		
		# Se não achou aresta, significa que é uma aresta reversa
		if not achouAresta:
			# Procura a aresta reversa
			for indiceAresta in range(G.shape[1]):
				if (G_aux[vertice][indiceAresta] == 1) and (G_aux[proximo][indiceAresta] == -1):
					arestas.append(indiceAresta * (-1)) # Adiciona a aresta com índice negativo
		

	return arestas



# Implementa de fato a visita DFS
def visitaDFS(G, visitados, s):
	
	# Marca o vértice atual como visitado
	visitados[s] = True
	for x in range(len(visitados)):
		if (not visitados[x]) and (G[s][x] > 0):
			visitaDFS(G, visitados, x)



# Função que encontra corte minimo
def corteMinimo(G, Gr):
	nArestas = G.shape[1] 	# Numero de arestas é o numero de colunas da matriz de incidência
	visitados = np.full((1, Gr.shape[1]), False).tolist()[0]

	# Procura o conjunto de todos os vertices alcançáveis a partir da origem (0)
	visitaDFS(Gr, visitados, 0)

	st_corte = np.zeros(nArestas)
	for aresta in range(nArestas):
		# Captura os indices dos vertices de chegada e saida da aresta
		u = int(np.where(G.getT()[aresta] == -1)[1])
		v = int(np.where(G.getT()[aresta] == 1)[1])

		''' Captura todas as arestas do grafo original que possuem valor 0 no residual
			Ou seja, todas as arestas que saem de vértices alcançáveis em Gr para vértices
			que não são alcançáveis (conjuntos U e V) '''
		if visitados[u] and not visitados[v]:
			if Gr[u][v] == 0:
				st_corte[aresta] = 1

	return st_corte



# ------------------------------------------------------------------------------------------- 
# PROGRAMA PRINCIPAL


# Recebe vetor c e matriz N da entrada
N = np.loadtxt(sys.argv[1], dtype = int, skiprows = 2)

# Extrai vetor c e matriz de incidência G como submatrizes
cT = np.asmatrix(N[0][:])
G = np.asmatrix(N[1:][:])


# TEMPORÁRIO >>> KIT PRINT DEBUG
#print("Custo das arestas: ", cT)
#print("\nMatriz de incidência (vertices x arestas): \n", G)


# Inicio do algoritmo de Ford-Fulkerson
# Inicializa o vetor com o fluxo de cada aresta (inicialmente é zero)
f = np.zeros_like(cT)
#print("\nFluxo inicial passando pelas arestas: ", f)


# Duplica matriz de incidência, capacidades e vetor de fluxos para inserir arestas reversas
#G_rev = np.asmatrix(N[0][:])
f_reversas = np.zeros_like(cT)

#G = np.concatenate(G_ind, G_rev)
#f = np.concatenate(f_inc, f_rev)

# Atualiza o grafo residual (definido como uma matriz de adjacência)
Gr = atualizaGrafoResidual(G, f, cT)
#print("\nGrafo residual inicial (matriz de adjacencia): \n", Gr)

# Definicao dos nós inicial e final para busca DFS (procura um caminho entre s-t)
s = 0	# Origem do fluxo é sempre o primeiro vértice
t = Gr.shape[0] - 1   # Destino do fluxo é sempre o último vértice


# Loop de execução do algoritmo
(continua, Cp, caminhoVertices) = achaCaminhoDFS(Gr, s, t) 	# Verifica se existe um caminho entre s-t
fluxoMaximo = 0
#print("\nCaminho encontrado? ", continua, ". Se sim, qual? ", caminhoVertices)

while continua:
	caminhoArestas = np.zeros(f.shape[1])
	arestas = achaArestasCaminho(G, caminhoVertices)
	
	# Adiciona fluxo encontrado às arestas	
	fluxoMaximo += Cp
	#print("Fluxo do caminho: ", Cp)
	for aresta in arestas:
		# Trata casos de arestas reversas
		if aresta < 0:
			aresta *= -1
			f.A[0][aresta] -= Cp
		# Se for arestas existentes no grafo original
		else:
			caminhoArestas[aresta] = 1
			f.A[0][aresta] += Cp
	
	#print("Caminho de arestas: ", arestas)

	# Impressão dos resultados parciais
	print(np.asarray(list(map(int, caminhoArestas.tolist()))))
	print(f.A[0])
	print(cT.A[0], "\n")

	#print("\nAtualizando informacoes...")
	#print("Fluxo atualizado por aresta: ", f.A[0])
	#print("Capacidade maxima das arestas: ", cT.A[0])

	Gr = atualizaGrafoResidual(G, f, cT)
	(continua, Cp, caminhoVertices) = achaCaminhoDFS(Gr, s, t)
	

	#print("\n\n-------------------------------------------------------------------")
	#print("\nNovo grafo residual (matriz de adjacencia): \n", Gr)
	#print("\nCaminho encontrado? ", continua, ". Se sim, qual? ", caminhoVertices)


# Impressão dos resultados finais
print("Fluxo maximo: ", int(fluxoMaximo))
print("Corte minimo: ", np.asarray(list(map(int, corteMinimo(G, Gr).tolist()))))
