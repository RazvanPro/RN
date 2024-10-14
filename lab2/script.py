import re
import math

def parse_system_from_file(filename):
    A = []
    B = []
    
    equation_pattern = re.compile(r'([-+]?\d*)x\s*([-+]\s*\d*)y\s*([-+]\s*\d*)z\s*=\s*([-+]?\d+)')

    with open(filename, 'r') as file:
        for line in file:
            match = equation_pattern.match(line.replace(" ", ""))
            if match:
                coeffs = [match.group(i) for i in range(1, 5)]

                row = [int(coeff if coeff not in ['', '+', '-'] else coeff + '1') for coeff in coeffs[:-1]]
                A.append(row)
                B.append(int(coeffs[-1])) 
    
    return A, B

def determinant_3x3(matrix):
    a11, a12, a13 = matrix[0]
    a21, a22, a23 = matrix[1]
    a31, a32, a33 = matrix[2]
    
    det = (a11 * (a22 * a33 - a23 * a32) 
           - a12 * (a21 * a33 - a23 * a31) 
           + a13 * (a21 * a32 - a22 * a31))
    return det

def trace(matrix):
    return matrix[0][0] + matrix[1][1] + matrix[2][2]

def vector_norm(vector):
    return math.sqrt(sum([x**2 for x in vector]))

def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def matrix_vector_multiplication(matrix, vector):
    return [sum(matrix[i][j] * vector[j] for j in range(len(vector))) for i in range(len(matrix))]

def cramer_solve(A, B):
    det_A = determinant_3x3(A)
    if det_A == 0:
        raise ValueError("Matrix A is singular, no unique solution exists.")

    def replace_column(matrix, col_index, vector):
        new_matrix = [row[:] for row in matrix] 
        for i in range(len(matrix)):
            new_matrix[i][col_index] = vector[i]
        return new_matrix

    Ax = replace_column(A, 0, B)
    Ay = replace_column(A, 1, B)
    Az = replace_column(A, 2, B)

    x = determinant_3x3(Ax) / det_A
    y = determinant_3x3(Ay) / det_A
    z = determinant_3x3(Az) / det_A

    return [x, y, z]

def cofactor(matrix, row, col):
    minor = [[matrix[i][j] for j in range(len(matrix)) if j != col] for i in range(len(matrix)) if i != row]
    return determinant_2x2(minor) if len(minor) == 2 else determinant_3x3(minor)

def determinant_2x2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

def adjugate(matrix):
    cofactor_matrix = [[((-1) ** (i + j)) * cofactor(matrix, i, j) for j in range(len(matrix))] for i in range(len(matrix))]
    return transpose(cofactor_matrix)

def inverse_matrix(matrix):
    det_A = determinant_3x3(matrix)
    if det_A == 0:
        raise ValueError("Matrix A is singular, no inverse exists.")
    
    adj_A = adjugate(matrix)
    return [[adj_A[i][j] / det_A for j in range(len(adj_A))] for i in range(len(adj_A))]

filename = 'equations.txt'
A, B = parse_system_from_file(filename)

print(f"Matrix A: {A}")
print(f"Vector B: {B}")

det_A = determinant_3x3(A)
print(f"Determinant of A: {det_A}")

trace_A = trace(A)
print(f"Trace of A: {trace_A}")

norm_B = vector_norm(B)
print(f"Norm of B: {norm_B}")

transpose_A = transpose(A)
print(f"Transpose of A: {transpose_A}")

AxB = matrix_vector_multiplication(A, B)
print(f"A * B = {AxB}")

solution_cramer = cramer_solve(A, B)
print(f"Solution using Cramer's rule: {solution_cramer}")

inv_A = inverse_matrix(A)
print(f"Inverse of A: {inv_A}")

solution_inversion = matrix_vector_multiplication(inv_A, B)
print(f"Solution using Matrix Inversion: {solution_inversion}")
