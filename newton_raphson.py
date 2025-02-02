
import math

def newton_raphson(f, df, p0, TOL, N=50):
    print("{:<10} {:<15} {:<15} ".format("Iteration", "po", "p1"))
    for i in range(N):
        if df(p0) == 0:
            print( "Derivative is zero at p0, method cannot continue.")
            return

        p = p0 - f(p0) / df(p0)

        if abs(p - p0) < TOL:
            return p
        print("{:<10} {:<15.9f} {:<15.9f} ".format(i, p0, p))
        p0 = p
    return p



if __name__ == "__main__":

    f = lambda x: (math.sin(x ** 2 + 5 * x + 6)) / (2 * math.exp(-x))
    df = lambda x: ((2 * x + 5) * math.cos(x ** 2 + 5 * x + 6)) / (2 * math.exp(-x)) + \
                   math.sin(x ** 2 + 5 * x + 6) / (2 * math.exp(-x))


    initial_guess = -2.5
    tol = 1e-6

    root_newton = newton_raphson(f, df, initial_guess, tol)
    print(f"Newton-Raphson Method: The root is approximately x = {root_newton:.6f}")
