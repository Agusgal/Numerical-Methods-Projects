import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Implementaciones.EqSist.SistemasEq import svd, leastSq




if __name__ == '__main__':

    df = pd.read_csv('p43.csv')
    x = np.array(df['x'].tolist())
    y = np.array(df['y'].tolist())
    n = x.shape[0]
    fig, ax = plt.subplots()
    ax.plot(x, y, '.', label='puntos')
    ax.legend()
    plt.show()

    A = np.zeros(shape=(len(x), 3))

    col1 = np.cos(np.sqrt(np.abs(x)))
    col2 = np.sin(np.sqrt(np.abs(x)))
    col3 = np.ones(len(x))

    A[:, 0] = col1
    A[:, 1] = col2
    A[:, 2] = col3

    sol = leastSq(A, y, method='qr')

    #U, S, V = svd(A)

    #xfinal = np.zeros(A.shape[1])

    #for k in range(A.shape[1]):
    #    xfinal += np.dot(U.T[k], y)/S[k][k] * V.T[k]

    def f(x):
        return sol[0] * np.cos(np.sqrt(abs(x))) + sol[1] * np.sin(np.sqrt(abs(x))) + sol[2]

    x2 = range(-50, 50)
    ig, ax = plt.subplots()
    ax.plot(x, y, '.', label='puntos')
    ax.plot(x2, [f(i) for i in x2], 'c', label="Función de ajuste")
    ax.legend()
    plt.show()

