import numpy as np

# === ROOT FINDING METHODS === #
def secant_method(f, x0, x1, tol=1e-6, max_iterations=50):
    for _ in range(max_iterations):
        f_x0 = f(x0)
        f_x1 = f(x1)
        denominator = f_x1 - f_x0

        if denominator == 0:
            return None  

        p = x0 - f_x0 * ((x1 - x0) / denominator)

        if abs(p - x1) < tol:
            return p 

        x0, x1 = x1, p

    return p

def newton_raphson(f, df, x0, tol=1e-6, max_iterations=50):
    for _ in range(max_iterations):
        derivative = df(x0)
        if derivative == 0:
            return None  

        x1 = x0 - f(x0) / derivative

        if abs(x1 - x0) < tol:
            return x1 

        x0 = x1

    return x0

# === INTEGRATION METHODS === #
def gaussian_quadrature(f, a, b, n=4):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    x = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    return np.sum(weights * f(x)) * 0.5 * (b - a)

def romberg_integration(f, a, b, max_iter=5, tol=1e-6):
    R = np.zeros((max_iter, max_iter))
    h = b - a
    R[0, 0] = 0.5 * h * (f(a) + f(b))

    for i in range(1, max_iter):
        h /= 2
        sum_term = np.sum([f(a + k * h) for k in range(1, 2 ** i, 2)])
        R[i, 0] = 0.5 * R[i - 1, 0] + h * sum_term

        for j in range(1, i + 1):
            R[i, j] = R[i, j - 1] + (R[i, j - 1] - R[i - 1, j - 1]) / ((4 ** j) - 1)

        if i > 1 and abs(R[i, i] - R[i - 1, i - 1]) < tol:
            return R[i, i]

    return R[max_iter - 1, max_iter - 1]

# === LINEAR SYSTEM SOLVERS === #
def gauss_seidel(A, b, x0, tol=1e-16, max_iter=200):
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    x = np.array(x0, dtype=float)

    for _ in range(max_iter):
        x_new = np.copy(x)
        for i in range(len(A)):
            sum1 = sum(A[i][j] * x_new[j] for j in range(i))
            sum2 = sum(A[i][j] * x[j] for j in range(i + 1, len(A)))
            x_new[i] = (b[i] - sum1 - sum2) / A[i][i]

        if np.linalg.norm(x_new - x, np.inf) < tol:
            return x_new
        x = x_new
    return x

def lu_factorization(A, b):
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    return np.linalg.solve(A, b)

def gauss_elimination(A, b):
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    return np.linalg.solve(A, b)

# === INTERPOLATION METHODS === #
def lagrange_interpolation(x_data, y_data, x_interpolate):
    n, result = len(x_data), 0

    for i in range(n):
        term = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x_interpolate - x_data[j]) / (x_data[i] - x_data[j])
        result += term

    return result

def neville_interpolation(x_data, y_data, x_interpolate):
    n = len(x_data)
    tableau = np.zeros((n, n))

    for i in range(n):
        tableau[i][0] = y_data[i]

    for j in range(1, n):
        for i in range(n - j):
            numerator = (x_interpolate - x_data[i+j]) * tableau[i][j-1] - (x_interpolate - x_data[i]) * tableau[i+1][j-1]
            denominator = x_data[i] - x_data[i+j]
            tableau[i][j] = numerator / denominator

    return tableau[0][n-1]

# === GENERAL SOLVER === #
def solve_problem(problem_number, problem_type, params):
    if problem_type == "root":
        return {
            "Problem": problem_number,
            "Secant Method": secant_method(params["f"], params["x0"], params["x1"], params["tol"], params["max_iter"]),
            "Newton-Raphson": newton_raphson(params["f"], params["df"], params["x0"], params["tol"], params["max_iter"])
        }
    elif problem_type == "integral":
        return {
            "Problem": problem_number,
            "Gaussian Quadrature": gaussian_quadrature(params["f"], params["a"], params["b"], 4),
            "Romberg Integration": romberg_integration(params["f"], params["a"], params["b"], params["max_iter"], params["tol"])
        }
    elif problem_type == "linear_system":
        return {
            "Problem": problem_number,
            "Gauss-Seidel": gauss_seidel(params["A"], params["b"], params["x0"], params["tol"], params["max_iter"]),
            "LU Factorization": lu_factorization(params["A"], params["b"])
        }
    elif problem_type == "linear_system_gauss":
        return {
            "Problem": problem_number,
            "Gauss-Seidel": gauss_seidel(params["A"], params["b"], params["x0"], params["tol"], params["max_iter"]),
            "Gauss Elimination": gauss_elimination(params["A"], params["b"])
        }
    elif problem_type == "interpolation":
        return {
            "Problem": problem_number,
            "Lagrange Interpolation": lagrange_interpolation(params["x_data"], params["y_data"], params["x_interpolate"]),
            "Neville Interpolation": neville_interpolation(params["x_data"], params["y_data"], params["x_interpolate"])
        }
    else:
        return None

# === PROBLEM DEFINITIONS === #
problems = [
    (1, "root", {
        "f": lambda x: (np.sin(x**2 + 5*x + 6)) / (2*np.exp(-x)),
        "df": lambda x: ((2*x + 5) * np.cos(x**2 + 5*x + 6)) / (2*np.exp(-x)),
        "x0": -3, "x1": 1, "tol": 1e-6, "max_iter": 100
    }),
    (20, "integral", {
        "f": lambda x: (x * np.exp(-x**2 + 5*x)) * (2*x**2 - 3*x - 5),
        "a": 0.5, "b": 1, "max_iter": 5, "tol": 1e-6
    }),
    (31, "linear_system", {
        "A": [[1, 0, -1], [-0.5, 1, -0.25], [1, -0.5, 1]],
        "b": [0.2, -1.425, 2], "x0": [0, 0, 0], "tol": 1e-6, "max_iter": 100
    }),
    (21, "linear_system_gauss", {
        "A": [[1, 1/2, 1/3], [1/2, 1/3, 1/4], [1/3, 1/4, 1/5]],
        "b": [1, 0, 0], "x0": [0, 0, 0], "tol": 1e-16, "max_iter": 200
    }),
    (37, "interpolation", {
        "x_data": [1.2, 1.3, 1.4, 1.5, 1.6],
        "y_data": [3.5095, 3.6984, 3.9043, 4.1293, 4.3756],
        "x_interpolate": 1.37
    }),
]

# Run each problem and print results
for problem in problems:
    problem_number, problem_type, params = problem
    results = solve_problem(problem_number, problem_type, params)
    print(f"Results for Problem {problem_number}: {results}")
