import numpy as np
from math import copysign, hypot


def leastsq(A, b, mode='Gram-Schmidt'):
    (rowsN, colsN) = np.shape(A)


    if mode == 'Gram-Schmidt':
        # Se inicializa matriz ortogonal vacia Q
        Q = np.empty([rowsN, rowsN])
        counter = 0

        #Computo de la matriz Q
        for a in A.T:
            u = np.copy(a)
            for i in range(0, counter):
                projection = np.dot(np.dot(Q[:, i].T, a), Q[:, i])
                u -= projection

            e = u / np.linalg.norm(u)
            Q[:, counter] = e

            counter += 1

        #Computo matriz triangular R
        R = np.dot(Q.T, A)

    elif mode == 'Householder':
        # Incializo ambas matrices, Q y R
        Q = np.identity(rowsN)
        R = np.copy(A)

        # Itero sobre cada subvector-columna y
        # computo matriz householder.
        for counter in range(rowsN - 1):
            x_ = R[counter:, counter]

            e = np.zeros_like(x_)
            e[0] = copysign(np.linalg.norm(x_), -A[counter, counter])
            u = x_ + e
            v = u / np.linalg.norm(u)

            Q_cnt = np.identity(rowsN)
            Q_cnt[counter:, counter:] -= 2.0 * np.outer(v, v)

            R = np.dot(Q_cnt, R)
            Q = np.dot(Q, Q_cnt.T)

    elif mode == 'Givens':
        #Inicializo ambas matrices
        Q = np.identity(rowsN)
        R = np.copy(A)

        #Itero sobre matriz triangular inferior
        (rows, cols) = np.tril_indices(rowsN, -1, colsN)
        for (row, col) in zip(rows, cols):

            # Computo matriz de Givens
            if R[row, col] != 0:
                r = hypot(R[col, col], R[row, col])
                c = R[col, col] / r
                s = -R[row, col] / r

                G = np.identity(rowsN)
                G[[col, row], [col, row]] = c
                G[row, col] = s
                G[col, row] = -s

                R = np.dot(G, R)
                Q = np.dot(Q, G.T)


    ##Tengo, Q y R, ahora computo Q1, R1

    Q1 = Q[:, :colsN]
    R1 = R[:colsN, :]

    #Resuelvo el sistema
    x = np.dot(np.dot(np.linalg.inv(R1), Q1.T), b)
    return (Q, R, x)


def test():
    A1 = np.array([[3, -6],
                  [4, -8],
                  [0, 1]], dtype=np.float64)

    b1 = np.array([[-1], [7], [2]], dtype=np.float64)

    A2 = np.array([[1, 1],
                   [4, 1],
                   [9, 1]], dtype=np.float64)

    b2 = np.array([[3.1],
                  [8.8],
                  [20.2]], dtype=np.float64)

    A3 = np.array([[3.02, -1.05, 2.53],
                   [4.33, 0.56, -1.78],
                   [-0.83, -0.54, 1.47]], dtype=np.float64)

    b3 = np.array([[-1.61], [7.23], [-3.38]]) ##Para este caso esta mal condicionada la matriz

    ##Con A1

    (QGS1, RGS1, testGS1) = leastsq(A1, b1)
    (QGiv1, RGiv1, testGivens1) = leastsq(A1, b1, 'Givens')
    (QHous1, RHous1, testH1) = leastsq(A1, b1, 'Householder')

    difGS1 = np.dot(QGS1, RGS1) - A1
    difGiv1 = np.dot(QGiv1, RGiv1) - A1
    difHous1 = np.dot(QHous1, RHous1) - A1

    print('x con matriz A1')
    print('Gram-Schmidt, x')
    prettyPrint(testH1)
    print('\n')

    print('Givens, x')
    prettyPrint(testGivens1)
    print('\n')

    print('Householder, x')
    prettyPrint(testH1)
    print('\n')

    print('QR - A, para matriz A1:\n')
    print('Gram-Schmidt')
    prettyPrint(difGS1)

    print('\n')

    print('Givens')
    prettyPrint(difGiv1)

    print('\n')
    print('Householder')
    prettyPrint(difHous1)
    print('\n')


    ##Caso con A2
    (QGS2, RGS2, testGS2) = leastsq(A2, b2)
    (QGiv2, RGiv2, testGivens2) = leastsq(A2, b2, 'Givens')
    (QHous2, RHous2, testH2) = leastsq(A2, b2, 'Householder')

    difGS2 = np.dot(QGS2, RGS2) - A2
    difGiv2 = np.dot(QGiv2, RGiv2) - A2
    difHous2 = np.dot(QHous2, RHous2) - A2

    print('x con matriz A2')
    print('Gram-Schmidt, x')
    prettyPrint(testGS2)
    print('\n')

    print('Givens, x')
    prettyPrint(testGS2)
    print('\n')

    print('Householder, x')
    prettyPrint(testH2)
    print('\n')


    print('QR - A, para matriz A2:\n')

    print('Gram-Schmidt')
    prettyPrint(difGS1)

    print('\n')
    print('Givens')
    prettyPrint(difGiv2)

    print('\n')
    print('HouseHolder')
    prettyPrint(difHous2)
    print('\n')

    ##Caso A3
    (QGS3, RGS3, testGS3) = leastsq(A3, b3)
    (QGiv3, RGiv3, testGivens3) = leastsq(A3, b3, 'Givens')
    (QHous3, RHous3, testH3) = leastsq(A3, b3, 'Householder')

    difGS3 = np.dot(QGS3, RGS3) - A3
    difGiv3 = np.dot(QGiv3, RGiv3) - A3
    difHous3 = np.dot(QHous3, RHous3) - A3

    print('x con matriz A3')
    print('Gram-Schmidt, x')
    prettyPrint(testGS3)
    print('\n')

    print('Givens, x')
    prettyPrint(testGS3)
    print('\n')

    print('Householder, x')
    prettyPrint(testH3)
    print('\n')

    print('QR - A, para matriz A3:\n')

    print('Gram-Schmidt')
    prettyPrint(difGS3)

    print('\n')
    print('Givens')
    prettyPrint(difGiv3)

    print('\n')
    print('HouseHolder')
    prettyPrint(difHous3)



#Funcion que imprime matrices mas bonitas UwU
def prettyPrint(mat):
    col_maxes = [max([len(("{:" + 'g' + "}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:" + str(col_maxes[i]) + 'g' + "}").format(y), end="  ")
        print("")

