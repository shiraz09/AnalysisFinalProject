import numpy as np

from colors import bcolors
from matrix_utility import swap_rows_elementary_matrix, row_addition_elementary_matrix


def lu(A):
    N = len(A)
    L = np.eye(N)

    for i in range(N):

        #
        pivot_row = i
        v_max = A[pivot_row][i]
        for j in range(i + 1, N):
            if abs(A[j][i]) > v_max:
                v_max = A[j][i]
                pivot_row = j


        if A[i][pivot_row] == 0:
            raise ValueError("can't perform LU Decomposition")


        if pivot_row != i:
            e_matrix = swap_rows_elementary_matrix(N, i, pivot_row)
            print(f"elementary matrix for swap between row {i} to row {pivot_row} :\n {e_matrix} \n")
            A = np.dot(e_matrix, A)
            print(f"The matrix after elementary operation :\n {A}")
            print(bcolors.OKGREEN, "---------------------------------------------------------------------------",
                  bcolors.ENDC)

        for j in range(i + 1, N):

            m = -A[j][i] / A[i][i]
            e_matrix = row_addition_elementary_matrix(N, j, i, m)
            e_inverse = np.linalg.inv(e_matrix)
            L = np.dot(L, e_inverse)
            A = np.dot(e_matrix, A)
            print(f"elementary matrix to zero the element in row {j} below the pivot in column {i} :\n {e_matrix} \n")
            print(f"The matrix after elementary operation :\n {A}")
            print(bcolors.OKGREEN, "---------------------------------------------------------------------------",
                  bcolors.ENDC)

    U = A
    return L, U



def backward_substitution(mat):
    N = len(mat)
    x = np.zeros(N)


    for i in range(N - 1, -1, -1):

        x[i] = mat[i][N]


        for j in range(i + 1, N):
            x[i] -= mat[i][j] * x[j]

        x[i] = (x[i] / mat[i][i])

    return x


def lu_solve(A_b):
    L, U = lu(A_b)
    print("Lower triangular matrix L:\n", L)
    print("Upper triangular matrix U:\n", U)

    result = backward_substitution(U)
    print(bcolors.OKBLUE, "\nSolution for the system:")
    for x in result:
        print("{:.6f}".format(x))


if __name__ == '__main__':

    A_b = [
        [1, 0, -1, 0.2],
        [-0.5, 1, -0.25, -1.425],
        [1, -0.5, 1, 2]
    ]

    print("Solving the system of linear equations using LU Decomposition:")
    lu_solve(A_b)
