from colors import bcolors
import numpy as np

def swap_row(mat, i, j):
    mat[i], mat[j] = mat[j], mat[i]

def gaussianElimination(mat):
    N = len(mat)

    singular_flag = forward_substitution(mat)

    if singular_flag != -1:

        if mat[singular_flag][N]:
            return "Singular Matrix (Inconsistent System)"
        else:
            return "Singular Matrix (May have infinitely many solutions)"

    return backward_substitution(mat)

def forward_substitution(mat):
    N = len(mat)
    for k in range(N):
        # Partial Pivoting
        pivot_row = k
        v_max = abs(mat[k][k])
        for i in range(k+1, N):
            if abs(mat[i][k]) > v_max:
                v_max = abs(mat[i][k])
                pivot_row = i

        if mat[pivot_row][k] == 0:
            return k


        if pivot_row != k:
            swap_row(mat, k, pivot_row)

        for i in range(k+1, N):
            m = mat[i][k] / mat[k][k]
            for j in range(k, N+1):
                mat[i][j] -= mat[k][j] * m
    return -1

def backward_substitution(mat):
    N = len(mat)
    x = np.zeros(N)


    for i in range(N-1, -1, -1):
        x[i] = mat[i][N]
        for j in range(i+1, N):
            x[i] -= mat[i][j] * x[j]
        x[i] /= mat[i][i]
    return x
if __name__ == '__main__':

    A_b = [
        [1,    1/2,  1/3,  1],
        [1/2,    1/3,    1/4,  0],
        [1/3,  1/4,  1/5,  0]
    ]


    mat_copy = [row[:] for row in A_b]


    result = gaussianElimination(mat_copy)

    if isinstance(result, str):

        print(bcolors.FAIL + result + bcolors.ENDC)
    else:

        print(bcolors.OKBLUE, "\nSolution for the system:")
        for idx, val in enumerate(result):
            print(f"x[{idx}] = {val:.6f}")
        print(bcolors.ENDC)
