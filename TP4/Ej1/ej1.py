import numpy as np
import matplotlib.pyplot as plt


def ruku4(f, t0, tf, deltaT, x0):

    # cantidad de puntos

    N = int((tf - t0)/deltaT)

    tk = np.linspace(t0, tf, N + 1)

    n = x0.shape[0]
    xk = np.zeros((N + 1, n))  # matriz del problema
    xk[0, :] = x0

    # Se itera sobre cada instante de tiempo
    for i in range(N):
        f1 = f(tk[i], xk[i, :])
        f2 = f(tk[i] + deltaT/2, xk[i, :] + ((deltaT/2) * f1))
        f3 = f(tk[i] + deltaT/2, xk[i, :] + ((deltaT/2) * f2))
        f4 = f(tk[i] + deltaT, xk[i, :] + deltaT * f3)

        xk[i + 1, :] = xk[i, :] + ((f1 + 2*f2 + 2*f3 + f4)/6) * deltaT

    return tk, xk





