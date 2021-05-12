import numpy as np
import matplotlib.pyplot as plt

def roku4(f, t0, tf, deltaT, x0):

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
        f4 = f(tk[i] + deltaT/2, xk[i, :] + deltaT * f3)

        xk[i + 1, :] = xk[i, :] + ((f1 + 2 * f2 + 2 * f3 + f4)/6) * deltaT

    return tk, xk

def dx(t, x):
    R = 1e3
    C = 1e-6
    W = 2.0 * np.pi * 1000
    A = 1.0
    T = 5*2*np.pi/W

    return ((A*np.cos(W*t)-x)/(R*C))


def known(time):
    R = 1e3
    C = 1e-6
    W = 2.0 * np.pi * 1000
    A = 1.0

    x = -np.exp(-time/(R*C))+np.cos(W*time)+W*R*C*np.sin(W*time)
    x = (A/(1+(W*R*C)**2))*x
    return x


if __name__ == '__main__':
    arr = np.array([[1, 2, 3], [4, 5, 6]])
    print(arr[0, :])






    x0 = np.zeros(1)
    W = 2.0 * np.pi * 1000
    T = 5 * 2 * np.pi / W
    h = T/10000

    t, xe = roku4(dx, 0, T, h, x0)
    xknown = known(t)

    fig, ax = plt.subplots()
    ax.plot(t, xknown, label='Conocida')
    ax.plot(t, xe[:, 0], label='Roku4')

    ax.legend()
    plt.show()


    fig, ax = plt.subplots()
    ax.plot(t, xe[:, 0] - xknown, label='Error')
    ax.legend()
    plt.show()


