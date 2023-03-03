import sympy
import sympy.physics.quantum.spin as ssympy
from sympy.physics.quantum import Operator, represent, Ket, qapply
from sympy import sqrt, Rational

import numpy as np

from timeit import default_timer as timer

from cefmatrices import *

def pretty_print(s):
    sympy.pretty_print(s, num_columns=220)  # For large screen in the office

def diagonalize(H):
    start_time = timer()

    P, D = H.diagonalize()

    end_time = timer()
    diag_str = f'Diagonalization time = {end_time-start_time} s'
    print(diag_str)

    return P, D

sympy.init_printing(use_unicode=False)
B20, B22, B40, B42, B44, B60, B62, B64, B66 = sympy.symbols('B20, B22, B40, B42, B44, B60, B62, B64, B66')
a11, a12, a13, a23, a22, a33 = sympy.symbols('a11, a12, a13, a23, a22, a33')
a,b,c,d,e = sympy.symbols('a,b,c,d,e')
m = sympy.symbols('m')


# Symbolic analysis of the crystal field effects in CeAuAl3

# Dont forget to edit these fields
Jval = Rational(5,2)
Jvalstr = f'{Jval:.1f}'.replace('.','p')
# filenameprefix = f'C:/Users/Stekiel/Documents/GitHub/mikibox/mikibox/crysfipy/symbolic/Symbolic-output/J_{Jvalstr}_B20_B22'

H = sympy.Matrix(   B20 * O_20(Jval) + B22 * O_22(Jval) + B40 * O_40(Jval) + B42 * O_42(Jval) + B44 * O_44(Jval) + B60 * O_60(Jval))
#H = sympy.Matrix( [[a11, a12, a13], [a12, a22, a23], [a13, a23, a33]] )

H = sympy.Matrix(
    [ [5*a + b, 0, d, 0, c, 0], 
    [0, -a - 3*b, 0, e, 0, c], 
    [d, 0, -4*a + 2*b, 0, e, 0], 
    [0, e, 0, -4*a + 2*b, 0, d], 
    [c, 0, e, 0, -a - 3*b, 0], 
    [0, c, 0, d, 0, 5*a + b]]
)

operations = [
    ('swap_rows',1,2,0),
    ('swap_rows',2,4,0),
    ('swap_rows',3,6,0),
    ('swap_cols',1,2,0),
    ('swap_cols',2,4,0),
    ('swap_cols',3,6,0),
] # J=2 case?
operations = [
    ('swap_rows',1,2,0),
    ('swap_rows',2,4,0),
    ('swap_cols',1,2,0),
    ('swap_cols',2,4,0),
    # ('add_rows',0,2,-Rational(1,3)*sympy.sqrt(5)),
    # ('add_rows',2,0,Rational(15,14)/sympy.sqrt(5)),
    # ('add_cols',2,0,-3/sympy.sqrt(5)),
    # ('add_cols',0,2,Rational(9,42)*sympy.sqrt(5))
] # J=2.5
operations=[]

diagonalize_hamiltonian = False
find_eigenvalues = True

print('Initial Hamiltonian')
pretty_print(H)
print(H)
if diagonalize_hamiltonian:
    P, D = diagonalize(H)
    print(D)
    for p in P:
        print(p)

if find_eigenvalues:
    eigenvalues = H.eigenvals()
    for eigenvalue, multiplicity in eigenvalues.items():
        print(f'Multiplicity {multiplicity}')
        print(f'eigenvalue   {eigenvalue.subs(e,0)}')

    eig1, eig2, eig3 = eigenvalues.keys()

    pretty_print(sympy.simplify((eig2-eig3).subs(e,0)))

for op,n1,n2,a in operations:
    if op=='swap_cols':
        H = H.elementary_col_op(op='n<->m', col1=n1, col2=n2)
    elif op=='swap_rows':
        H = H.elementary_row_op(op='n<->m', row1=n1, row2=n2)
    elif op=='add_cols':
        H = H.elementary_col_op(op='n->n+km', col1=n1, col2=n2, k=a)
    elif op=='add_rows':
        H = H.elementary_row_op(op='n->n+km', row1=n1, row2=n2, k=a)
    else:
        raise ValueError('Invalid matrix operation')

    print('New Hamiltonian')
    H = sympy.simplify(H)
    pretty_print(H)

H1 = H[0:4,0:4]
# pretty_print(H1, wrap_line=True)

if False:
    diagonalize(H1)
#print(H.eigenvals())
