import numpy as np


def hhDer(time, x):
    v0 = 0.50

    s = x[0]
    p = x[1]

    ds = v0 - 0.23*s * p**2
    dp = 0.23*s * p**2 - 0.4*p

    dx = np.zeros_like(x)
    dx[0] = ds
    dx[1] = dp

    return dx

