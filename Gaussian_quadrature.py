import math
from math import sqrt
import numpy as np

from colors import bcolors


def gaussian_quadrature(f, a, b, n):
    """
    Gaussian Quadrature for Numerical Integration

    Parameters:
    f (function): The function to be integrated.
    a (float): The lower limit of integration.
    b (float): The upper limit of integration.
    n (int): The number of nodes and weights to use.

    Returns:
    float: The approximate definite integral of the function over [a, b].
    """

    nodes, weights = np.polynomial.legendre.leggauss(n)


    x = 0.5 * (b - a) * nodes + 0.5 * (a + b)


    integral = sum(weights[i] * f(x[i]) for i in range(n))
    integral *= 0.5 * (b - a)

    return integral


def sub_division(f,a,b,tol,entire,results):

    a_z = a+(b-a)/2.
    b_k = a+(b-a)/2.
    entire=gaussian(f,a,b)
    right=gaussian(f,a_z,b)
    left=gaussian(f,a,b_k)
    if abs(entire-(left+right))<tol * max(abs(entire), (abs(left)+abs(right))):
        results.append(entire)
        return entire
    x=sub_division(f,a_z,b,tol,right,results)+sub_division(f,a,b_k,tol,left,results)
    results.append(x)
    return x

def gaussian(f,a,b):
    u=(b-a)/2.*(5./9*f((b-a)/2.*-1.*sqrt(3./5)+(b+a)/2.)+8./9*f((b+a)/2.)+5./9*f((b-a)/2.*sqrt(3./5)+(b+a)/2.))
    return u

def adaptive_gaussian_quadrature(f,a,b,tol,results):

    return sub_division(f,a,b,tol,gaussian(f,a,b),results)


if __name__ == '__main__':

    f = lambda x: (x * math.exp(-x ** 2 + 5 * x)) * (2 * x ** 2 - 3 * x - 5)
    a = 0.5
    b = 1
    tol = 1e-6

    results = []


    res = adaptive_gaussian_quadrature(f, a, b, tol, results)

    for i, result in enumerate(results):
        print(f"- Result number {i + 1}: {result}")


    print(bcolors.OKBLUE, f"\nIntegral value in range [{a},{b}] is {res}", bcolors.ENDC)