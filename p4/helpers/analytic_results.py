import numpy as np
import sympy as sym
from sympy import *
from IPython import embed


if __name__=="__main__":
    sigma_tt = lambda r: 1 / 3 * (1 + 4 / r ** 2)
    sigma_rr = lambda r: 1 / 3 * (1 - 4 / r ** 2)

    rs = np.linspace(1.05, 1.95, 10)

    for r in rs:
        print(f"r = {r}")
        print(f"sigma_tt = {sigma_tt(r)}")
        print(f"sigma_rr = {sigma_rr(r)}")
        print("\n")

    z = symbols('z')
    eq1 = (0.1 - z) / 0.1
    eq2 = z / 0.1

    int1 = integrate(eq1, (z, 0, 0.1))
    int2 = integrate(eq2, (z, 0, 0.1))

    embed()
