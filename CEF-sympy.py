import sympy
import numpy as np
from timeit import default_timer as timer
from cefmatrices import *

# Some wrapper funcitons
def pretty_print(s):
    sympy.pretty_print(s, num_columns=220)  # For large screen in the office, modify accordgin to the size of monitor

def diagonalize(H):
    start_time = timer()

    P, D = H.diagonalize()

    end_time = timer()
    diag_str = f'Diagonalization time = {end_time-start_time} s'
    print(diag_str)

    return P, D

# Initialization of stuff for sympy
sympy.init_printing(use_unicode=False)
B20, B22, B40, B42, B44, B60, B62, B64, B66 = sympy.symbols('B20, B22, B40, B42, B44, B60, B62, B64, B66')

### Main controls for the calculations
# Main variables
diagonalize_hamiltonian = False
find_eigenvalues = True
Jval = sympy.Rational(5,2)  # This can be any positive integer or half integer, like 1/2, 1, 3/2, 2 ...

### DEFINE THE HAMILTONIAN
#H = sympy.Matrix(   B20 * O_20(Jval) + B22 * O_22(Jval) + B40 * O_40(Jval) + B42 * O_42(Jval) + B44 * O_44(Jval) + B60 * O_60(Jval))
H = sympy.Matrix(   B20 * O_20(Jval) + B22 * O_22(Jval) + B40 * O_40(Jval))

print('Initial Hamiltonian')
pretty_print(H)
# print(H)

if diagonalize_hamiltonian:
    P, D = diagonalize(H)
    print(D)
    for p in P:
        print(p)

if find_eigenvalues:
    eigenvalues = H.eigenvals()
    for eigenvalue, multiplicity in eigenvalues.items():
        print(f'Multiplicity {multiplicity}')
        print(f'eigenvalue   {eigenvalue}')

    eig1, eig2, eig3 = eigenvalues.keys()

    pretty_print(sympy.simplify((eig2-eig3)))