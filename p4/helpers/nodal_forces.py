import numpy as np
import sympy as sym
from sympy import *
from IPython import embed


if __name__=="__main__":
    x, x1, x2 = symbols('x x1 x2')

    eq1 = (x - x2) / (x1 - x2)
    eq2 = (x - x1) / (x2 - x1)

    integral1 = sym.simplify(integrate(eq1, (x, x1, x2)))
    integral2 = sym.simplify(integrate(eq2, (x, x1, x2)))

    embed()
