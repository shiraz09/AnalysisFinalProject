
import math


def secant_method(f, x0, x1, TOL, N=50):
    print("{:<10} {:<15} {:<15} {:<15}".format("Iteration", "xo", "x1", "p"))
    for i in range(N):
        if f(x1) - f(x0) == 0:
            print( " method cannot continue.")
            return

        p = x0 - f(x0) * ((x1 - x0) / (f(x1) - f(x0)))

        if abs(p - x1) < TOL:
            return p  
        print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f}".format(i, x0, x1,p))
        x0 = x1
        x1 = p
    return p


if __name__ == "__main__":

    f = lambda x: math.sin(x**2 + 5*x + 6) / (2 * math.exp(-x))


    x0, x1 = -3, -2.5
    tol = 1e-6


    root_secant = secant_method(f, x0, x1, tol)
    if root_secant is not None:
        print(f"Secant Method: The root is approximately x = {root_secant:.6f}")
    else:
        print("Secant Method failed to converge.")
