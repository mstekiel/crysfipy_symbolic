import sympy
from sympy import pprint
import numpy as np

from cef_hamiltonian import CEF_Hamiltonian, Hamiltonian
# import mikibox as ms

from timeit import default_timer as timer


# Dont forget to edit these fields
# Jval = 2.5
# lattice='tetragonal_2'
# Ex.1 interesting set of 4D matrices
# Jval = 12
# lattice='hexagonal_2'
# Ex.2 5D matrix
# Jval = 12
# lattice='hexagonal_2'

# calculate_eigenvectors = False

# diagonalizeH = False
# investigateSubspaces = False

# exportLatex = False
# exportSympy = False

hline = '--------------------------------------------------------------------------------------'

Jval = 3
lattice='orthorhombic'
determine_eigenvectors = False

H = CEF_Hamiltonian(symmetry=lattice, Jval=Jval)

start_time = timer()

H.make_block_form()

H.determine_eigenvalues()
eigenvalues_timer = timer()
print(f'# Determining eigenvalues = {eigenvalues_timer-start_time:4f} s')

if determine_eigenvectors:
    H.determine_eigenvectors()
    eigenvectors_timer = timer()
    print(f'# Determining eigenvectors = {eigenvectors_timer-start_time:4f} s')



# # Evaluations
# sqrt  = np.emath.sqrt
# I = 1j

# CePdAl3
# Bij_names = ['B20', 'B22', 'B40', 'B42', 'B44', 'B60', 'B62', 'B64', 'B66']
# Bij_values = [1.203, -0.3   , -0.001, -1.1182781e-02,  0.244, 1, 0.01, 0.3, -0.001]    # orth

# Bij = ms.crysfipy.CEFpars('C4', [1.203, -0.001, 0.244], 'meV')
# Ce_4mm = ms.crysfipy.CEFion(ms.crysfipy.Ion('Ce'), (0,0,0), Bij)

# print(Ce_4mm)

# PrOs4Sb12
# doi.org/10.1103/PhysRevLett.93.157003
# Bij_names = ['B40','B60','B66']
# Bij_values = [0.2e-2, 0.11e-3, -0.9e-3]    # orth

# Bij_subs = {key:val for key, val in zip(Bij_names, Bij_values)}

# def cal_symbolic(Elevels, Bij_names, Bij_values):
#     e_symbolic = []

#     for E in Elevels:
#         E = str(E)
#         for Bij_name, Bij_value in zip(Bij_names, Bij_values):
#             E = E.replace(Bij_name, f'({Bij_value})')

#         e = eval(E)
#         if np.abs(np.imag(e))>1e10:
#             e = np.real(e)

#         e_symbolic.append(e)

#     return np.array(e_symbolic)-min(e_symbolic)

###
# Calculations

# H = CEF_Hamiltonian(symmetry=lattice, Jval=Jval)

# H.save_latex('./Symbolic-output/test.tex')

# print(H)
# H.make_block_form()
# for Hsub in H.subs:
#     print(Hsub)
#     # print(Hsub.print_latex())

# print('Main eigenvalues')
# print(H.determine_eigenvalues())
# print(H.determine_eigenvectors())

# print(hline)

# for Hsub in H.subs:
#     print(Hsub)
#     print(Hsub.eigenvalues)
#     print(Hsub.eigenvectors)

    # if Hsub.matrix.shape[0] < 4:
    #     print('Look into eigenvectors of this Hsub')
    #     Hsub.determine_eigenvectors()
    #     print(Hsub.eigenvectors)

# print(H.Neval_eigenvalues(Bij_subs, precision=5))


# print(hline)
# for Hsub in H.make_block_form():
#     print(Hsub)
#     eigenvalues = Hsub.determine_eigenvalues()
#     eigenvectors = Hsub.determine_eigenvectors()

#     print(eigenvalues)
#     print(eigenvectors)

#     print(Hsub.Neval_eigenvectors(Bij_subs, normalize=True))
#     evals_reformat = [formula for formula, _ in eigenvalues]
#     print(cal_symbolic(evals_reformat, Bij_names, Bij_values))
        
# N = 100
# print('Time numerical')
# test = timeit.Timer(lambda: Hsub.Neval_eigenvalues(Bij_subs)) 
# print (test.timeit(N))

# print(hline)
#################################
# Test Hamiltonian
# h = Hamiltonian(base=['1', '0', '-1'], matrix=sympy.eye(3))

# print(h)

# eigenvalues = h.determine_eigenvalues()
# eigenvectors = h.determine_eigenvectors()

# print(eigenvalues)
# print(eigenvectors)
