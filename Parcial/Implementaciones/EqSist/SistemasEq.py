import numpy as np
from math import sqrt
from scipy import linalg


# Eliminacion Gaussiana
def gaussian_elimination(A, b):
    U = np.copy(A)
    newb = np.copy(b)

    print(f'Inicio:\n{U}')  # comentar si no se quiere ver cada paso de la eliminacion
    # k: numero de paso de la eliminacion
    for k in range(len(U)):
        # j: fila que se esta eliminando
        for j in range(k + 1, len(U)):
            m = U[j][k] / U[k][k]
            U[j] = U[j] - m * U[k]
            newb[j] = newb[j] - m * newb[k]
        print(f'Paso {k + 1}:\n{U}')  # comentar si no se quiere ver cada paso de la eliminacion
    return (U, newb)


# Descomposición LU
def lu_decomposition(A):
    """(I | A) -> (L | A)"""
    U = np.copy(A)
    L = np.eye(len(A))
    print(f'Inicio:\nL:\n{L}\nU:\n{U}')  # comentar si no se quiere ver cada paso de la descomposicion
    # k: numero de paso de la descomposicion
    for k in range(len(U)):
        # j: fila que se esta eliminando
        for j in range(k+1, len(U)):
            m = U[j][k]/U[k][k]
            U[j] = U[j] - m*U[k]
            L[j][k] = m
        print(f'Paso {k + 1}:\nL:\n{L}\nU:\n{U}') # comentar si no se quiere ver cada paso de la descomposicion
    return (L, U)


#Factorización de Cholesky
def cholesky(A):
    L = [[0.0] * len(A) for _ in range(len(A))]
    for i, (Ai, Li) in enumerate(zip(A, L)):
        for j, Lj in enumerate(L[:i + 1]):
            s = sum(Li[k] * Lj[k] for k in range(j))
            Li[j] = sqrt(Ai[i] - s) if (i == j) else \
                (1.0 / Lj[j] * (Ai[j] - s))
    return L

#Decomposición QR
def qr(A):
    (rowsN, colsN) = np.shape(A)

    # Se inicializa matriz ortogonal vacia Q
    Q = np.empty([rowsN, rowsN])
    R = np.empty([rowsN, colsN])
    counter = 0

    # Computo de la matriz Q
    for a in A.T:
        u = np.copy(a)
        for i in range(0, counter):
            projection = np.dot(np.dot(Q[:, i].T, a), Q[:, i])
            u -= projection

        e = u / np.linalg.norm(u)
        Q[:, counter] = e

        counter += 1

    R = np.matmul(Q.T, A)

    Q1 = Q[:, :colsN]
    R1 = R[:colsN, :colsN]


    return Q1, R1


# Sustitución hacia atrás ----> A.x = b con A triangular superior
def solveBackSub(A, b):
    sol = np.zeros(shape=b.shape)
    n = len(b)

    sol[-1] = b[-1] / A[n - 1][n - 1]
    for k in range(n - 2, -1, -1):
        sol[k] = (b[k] - np.dot(A[k, k+1:n], sol[k+1:n]))/A[k, k]

    return sol


# Sustitucion hacia adelante ----> A.x = b con A triangular inferior
def solveForwSub(A, b):
    sol = np.zeros(shape=b.shape)
    n = len(b)

    sol[0] = b[0] / A[0][0]
    for k in range(1, n):
        sol[k] = (b[k] - np.dot(A[k, :k], sol[:k]))/A[k][k]

    return sol

# Descomposicion en valores singulares
def svd(A):

    (h, w), B, mat_range = A.shape, A.T.dot(A), min(A.shape)
    eigval, V = np.linalg.eig(B)

    idx = eigval.argsort()[::-1]
    eigval, V = np.sqrt(np.abs(eigval[idx])), np.real(V[:, idx])


    ##Calculo U a partir de V, Ui =1/eigVal * A * Vi
    U = np.zeros((h, h))
    U[:h, : mat_range] = np.dot(A, V[:, :mat_range]) / eigval[:mat_range]



    sigma = np.zeros(A.shape)
    np.fill_diagonal(sigma, eigval[:mat_range])

    return U, sigma, V


# Least squares
def leastSq(A, y, method='svd'):
    if method == 'svd':
        U, S, V = svd(A)

        xfinal = np.zeros(A.shape[1])

        for k in range(A.shape[1]):
            xfinal += np.dot(U.T[k], y) / S[k][k] * V.T[k]

        return xfinal

    elif method == 'qr':## ojo con data type de array
        Q1, R1 = qr(A)
        xfinal = np.linalg.solve(R1, np.matmul(Q1.T, y))
        return xfinal
    elif method == 'normalEqs':
        pass


if __name__ == '__main__':
    pass
