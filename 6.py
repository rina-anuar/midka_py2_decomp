import tkinter as tk
from tkinter import messagebox
import statistics

def matrix_multiply(A, B):
    return [[round(sum(a * b for a, b in zip(A_row, B_col)), 1)
             for B_col in zip(*B)] for A_row in A]

def transpose(matrix):
    return list(map(list, zip(*matrix)))

def identity_matrix(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def jacobi_rotation(A, p, q):
    n = len(A)
    R = identity_matrix(n)

    if abs(A[p][q]) < 1e-10:
        return R

    theta = (A[q][q] - A[p][p]) / (2 * A[p][q])
    t = 1 / (abs(theta) + (1 + theta**2)**0.5)
    if theta < 0:
        t = -t

    c = round(1 / (1 + t**2)**0.5, 1)
    s = round(t * c, 1)

    R[p][p] = R[q][q] = c
    R[p][q] = s
    R[q][p] = -s

    return R

def svd_decomposition(A, max_iterations=100, tolerance=1e-10):
    m, n = len(A), len(A[0])
    U = identity_matrix(m)
    V = identity_matrix(n)

    ATA = matrix_multiply(transpose(A), A)

    for _ in range(max_iterations):
        max_element = 0
        p, q = 0, 1

        for i in range(n):
            for j in range(i + 1, n):
                if abs(ATA[i][j]) > max_element:
                    max_element = abs(ATA[i][j])
                    p, q = i, j

        if max_element < tolerance:
            break

        J = jacobi_rotation(ATA, p, q)
        V = matrix_multiply(V, J)
        ATA = matrix_multiply(matrix_multiply(transpose(J), ATA), J)

    S = [0] * n
    for i in range(n):
        S[i] = round((abs(ATA[i][i]))**0.5, 1)

    for j in range(n):
        if S[j] > tolerance:
            for i in range(m):
                U[i][j] = round(sum(A[i][k] * V[k][j] for k in range(n)) / S[j], 1)

    # Normalize columns of U
    for j in range(n):
        norm = round((sum(U[i][j]**2 for i in range(m)))**0.5, 1)
        if norm > tolerance:
            for i in range(m):
                U[i][j] = round(U[i][j] / norm, 1)

    return U, S, transpose(V)

def lu_decomposition(A):
    n = len(A)
    L = [[0.0 for _ in range(n)] for _ in range(n)]
    U = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for k in range(i, n):
            sum_u = sum(L[i][j] * U[j][k] for j in range(i))
            U[i][k] = round(A[i][k] - sum_u, 1)

        L[i][i] = 1.0
        for k in range(i + 1, n):
            sum_l = sum(L[k][j] * U[j][i] for j in range(i))
            if U[i][i] != 0:
                L[k][i] = round((A[k][i] - sum_l) / U[i][i], 1)

    return L, U

def on_calculate():
    try:
        size = int(size_var.get())

        A = []
        for i in range(size):
            row = list(map(float, matrix_entries[i].get().split(',')))
            if len(row) != size:
                raise ValueError("The number of elements in the row must match the size of the matrix.")
            A.append(row)

        decomp_type = decomp_var.get()
        if decomp_type == "SVD":
            U, S, VT = svd_decomposition(A)
            result = (U, S, VT)
            U_str = '\n'.join([' '.join(map(str, row)) for row in U])
            S_str = ' '.join(map(str, S))
            VT_str = '\n'.join([' '.join(map(str, row)) for row in VT])

            # Calculate the standard deviation of singular values
            std_deviation = round(statistics.stdev(S), 2) if len(S) > 1 else 0.0

            S_matrix = [[S[i] if i == j else 0 for j in range(len(S))] for i in range(len(S))]
            reconstructed = matrix_multiply(matrix_multiply(U, S_matrix), VT)
            reconstructed_str = '\n'.join([' '.join(map(str, row)) for row in reconstructed])

            output_label.config(text=f"SVD decomposition:\nU:\n{U_str}\n\nSingular values:\n{S_str}\n\nV^T:\n{VT_str}\n\nStandard deviation of singular values: {std_deviation}\n\nReconstructed matrix:\n{reconstructed_str}")
        elif decomp_type == "LU":
            L, U = lu_decomposition(A)
            result = (L, U)
            output_label.config(text=f"LU decomposition:\nL: {L}\nU: {U}")
    except Exception as e:
        messagebox.showerror("Error", str(e))

def on_size_change(*args):
    size = int(size_var.get())
    for entry in matrix_entries:
        entry.grid_remove()
    for i in range(size):
        matrix_entries[i].grid(row=i + 2, column=0)

root = tk.Tk()
root.title("Matrix Decomposition")
size_var = tk.StringVar(value="2")
size_label = tk.Label(root, text="Select matrix size (2x2 to 5x5):")
size_label.grid(row=0, column=0)

size_options = [2, 3, 4, 5]
size_menu = tk.OptionMenu(root, size_var, *size_options, command=on_size_change)
size_menu.grid(row=0, column=1)

matrix_entries = [tk.Entry(root) for _ in range(5)]
for i in range(5):
    matrix_entries[i].grid(row=i + 2, column=0)

decomp_var = tk.StringVar(value="SVD")
decomp_label = tk.Label(root, text="Select decomposition:")
decomp_label.grid(row=7, column=0)

decomp_menu = tk.OptionMenu(root, decomp_var, "SVD", "LU")
decomp_menu.grid(row=7, column=1)

calculate_button = tk.Button(root, text="Calculate", command=on_calculate)
calculate_button.grid(row=8, column=0, columnspan=2)

output_label = tk.Label(root, text="")
output_label.grid(row=9, column=0, columnspan=2)

on_size_change()

root.mainloop()