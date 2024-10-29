import tkinter as tk
from tkinter import messagebox

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

def qr_decomposition(A):
    n = len(A)
    m = len(A[0])
    Q = [[A[i][j] for j in range(m)] for i in range(n)]
    R = [[0.0] * m for _ in range(m)]
    
    for j in range(m):
        R[j][j] = round((sum(Q[i][j]**2 for i in range(n)))**0.5, 1)
        
        for i in range(n):
            Q[i][j] = round(Q[i][j] / R[j][j], 1) if R[j][j] != 0 else 0
        
        for k in range(j + 1, m):
            R[j][k] = round(sum(Q[i][j] * A[i][k] for i in range(n)), 1)
            for i in range(n):
                Q[i][k] = round(Q[i][k] - Q[i][j] * R[j][k], 1)
    
    return Q, R

def verify_decomposition(A, decomp_type, matrices):
    if decomp_type == "SVD":
        U, S, VT = matrices
        S_matrix = [[S[i] if i == j else 0 for j in range(len(S))] for i in range(len(S))]
        reconstructed = matrix_multiply(matrix_multiply(U, S_matrix), VT)
    elif decomp_type == "LU":
        L, U = matrices
        reconstructed = matrix_multiply(L, U)
    else:  # QR
        Q, R = matrices
        reconstructed = matrix_multiply(Q, R)
    
    print(f"\nChecking {decomp_type}-decomposition:")
    print("Restored matrix:")
    for row in reconstructed:
        print([round(x, 1) for x in row])
    
    print("\nOriginal matrix:")
    for row in A:
        print([round(x, 1) for x in row])

def print_matrices(name, matrices):
    print(f"\n{name} decomposition:")
    for label, matrix in matrices:
        print(f"\n{label}:")
        if isinstance(matrix, list) and isinstance(matrix[0], list): 
            for row in matrix:
                print([round(x, 1) for x in row])
        else:  
            print([round(x, 1) for x in matrix])

def on_calculate():
    try:
        size = int(size_var.get())
        
        A = []
        for i in range(size):
            row = list(map(float, matrix_entries[i].get().split(',')))
            if len(row) != size:
                raise ValueError("The number of elements in a row must match the size of the matrix.")
            A.append(row)
        
        decomp_type = decomp_var.get()
        if decomp_type == "SVD":
            U, S, VT = svd_decomposition(A)
            result = (U, S, VT)
            U_str = '\n'.join([' '.join(map(str, row)) for row in U])
            S_str = ' '.join(map(str, S))
            VT_str = '\n'.join([' '.join(map(str, row)) for row in VT])
            
            S_matrix = [[S[i] if i == j else 0 for j in range(len(S))] for i in range(len(S))]
            reconstructed = matrix_multiply(matrix_multiply(U, S_matrix), VT)
            reconstructed_str = '\n'.join([' '.join(map(str, row)) for row in reconstructed])
            
            output_label.config(text=f"SVD decomposition:\nU:\n{U_str}\n\nSingular values:\n{S_str}\n\nV^T:\n{VT_str}\n\nRestored matrix:\n{reconstructed_str}")
            verify_decomposition(A, "SVD", result)
        elif decomp_type == "LU":
            L, U = lu_decomposition(A)
            result = (L, U)
            L_str = '\n'.join([' '.join(map(str, row)) for row in L])
            U_str = '\n'.join([' '.join(map(str, row)) for row in U])
            output_label.config(text=f"LU decomposition:\nL:\n{L_str}\n\nU:\n{U_str}")
            verify_decomposition(A, "LU", result)
        elif decomp_type == "QR":
            Q, R = qr_decomposition(A)
            result = (Q, R)
            Q_str = '\n'.join([' '.join(map(str, row)) for row in Q])
            R_str = '\n'.join([' '.join(map(str, row)) for row in R])
            output_label.config(text=f"QR decomposition:\nQ:\n{Q_str}\n\nR:\n{R_str}")
            verify_decomposition(A, "QR", result)
    except Exception as e:
        messagebox.showerror("Ошибка", str(e))

def on_size_change(*args):
    size = int(size_var.get())
    for entry in matrix_entries:
        entry.grid_remove()
    for i in range(size):
        matrix_entries[i].grid(row=i + 2, column=0)

root = tk.Tk()
root.title("Decomposition matrix")

size_var = tk.StringVar(value="2")
size_label = tk.Label(root, text="Select matrix size (from 2x2 to 5x5):")
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

decomp_menu = tk.OptionMenu(root, decomp_var, "SVD", "LU", "QR")
decomp_menu.grid(row=7, column=1)

calculate_button = tk.Button(root, text="Calculate", command=on_calculate)
calculate_button.grid(row=8, column=0, columnspan=2)

output_label = tk.Label(root, text="")
output_label.grid(row=9, column=0, columnspan=2)

on_size_change()

root.mainloop()