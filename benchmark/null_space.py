from random import seed, choice
import numpy as np
import galois
import re


seed(4269)


N = 1000
M = [[choice([0,1]) for i in range(N)] for j in range(N+1)]

def original(M):
	def _to_rref(M_):
		n,m = M_.shape
		marked = []
		for j in range(m):
			i = 0
			while i < n and M_[i,j] != 1:
				i += 1
			if i == n:
				continue
			else:
				marked.append(i)
				for k in range(m):
					if k != j:
						if M_[i, k] == 1:
							M_[:, k] ^= M_[:, j]
		return marked
	marked = _to_rref(M)
	dependent_index = [i for i in range(len(M)) if i not in marked][0]
	basis = [dependent_index]
	ones_cols = np.where(M[dependent_index] == 1)[0]
	basis.extend([ np.where(M[:, j] == 1)[0][0] for j in ones_cols ])
	return basis

def w_galois(M):
    gf = galois.GF(2)

# print(original(np.array(M)))
# print(M)
