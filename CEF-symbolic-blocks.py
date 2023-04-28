import sympy
import numpy as np

from timeit import default_timer as timer

from cefmatrices import *

# sympy.init_printing(use_unicode=False)
B20, B22, B2m2 = sympy.symbols('B20, B22, B2m2', real=True)
B40, B42, B4m2, B43, B4m3, B44, B4m4 = sympy.symbols('B40, B42, B4m2, B43, B4m3, B44, B4m4', real=True)
B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6 = sympy.symbols('B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6', real=True)


# Dont forget to edit these fields
Jval = 2.5
lattice='orthorhombic'
# Ex.1 interesting set of 4D matrices
# Jval = 11.5   
# lattice='hexagonal_2'
# Ex.2 5D matrix
# Jval = 12
# lattice='hexagonal_2'

diagonalizeH = False
investigateSubspaces = True

exportLatex = True
exportSympy = False

###
# Calculations

Jvalstr = f'{Jval:.1f}'.replace('.','p')
filenameprefix = f'./Symbolic-output/J_{Jvalstr}_{lattice}'


# Prepare a list of kets that form the basis
if Jval%1==0:
    freeionkets = [f'|{int(x)}>' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]
elif Jval%1==0.5:
    freeionkets = [f'|{int(2*x):d}/2>' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]

# Symmetry allowed Bnm according to McPhase manual
# https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html
H = sympy.Matrix({
'monoclinic':
    B20*O_20(Jval) + B22*O_22(Jval) + B2m2*O_2m2(Jval) + \
    B40*O_40(Jval) + B42*O_42(Jval) + B4m2*O_4m2(Jval) + B44*O_44(Jval) + B4m4*O_4m4(Jval) + \
    B60*O_60(Jval) + B62*O_62(Jval) + B6m2*O_6m2(Jval) + \
    B64*O_64(Jval) + B6m4*O_6m4(Jval) + B66*O_66(Jval) + B6m6*O_6m6(Jval),
'orthorhombic':
    B20*O_20(Jval) + B22*O_22(Jval) + \
    B40*O_40(Jval) + B42*O_42(Jval) + B44*O_44(Jval) + \
    B60*O_60(Jval) + B62*O_62(Jval) + B64*O_64(Jval) + B66*O_66(Jval),
'tetragonal_1':
    B20*O_20(Jval) + \
    B40*O_40(Jval) + B44*O_44(Jval) + B4m4*O_4m4(Jval) + \
    B60*O_60(Jval) + B64*O_64(Jval) + B6m4*O_6m4(Jval),
'tetragonal_2':
    B20*O_20(Jval) + B40*O_40(Jval) + B44*O_44(Jval) + \
    B60*O_60(Jval) + B64*O_64(Jval),
'trigonal_1':
    B20*O_20(Jval) + \
    B40*O_40(Jval) + B43*O_43(Jval) + B4m3*O_4m3(Jval) + \
    B60*O_60(Jval) + B63*O_63(Jval) + B6m3*O_6m3(Jval) + B66*O_66(Jval) + B6m6*O_6m6(Jval),
'trigonal_2':
    B20*O_20(Jval) + B40*O_40(Jval) + B43*O_43(Jval) + \
    B60*O_60(Jval) + B63*O_63(Jval) + B66*O_66(Jval),
'hexagonal_1':
    B20*O_20(Jval) + B40*O_40(Jval) + \
    B60*O_60(Jval) + B66*O_66(Jval) + B6m6*O_6m6(Jval),
'hexagonal_2':
    B20*O_20(Jval) + B40*O_40(Jval) + B60*O_60(Jval) + B66*O_66(Jval),
'cubic_1':
    B40*O_40(Jval) + sympy.Rational(5,2)*B40*O_44(Jval) + \
    B60*O_60(Jval) + -1*B66*O_62(Jval) + sympy.Rational(-21, 2)*B60*O_64(Jval) + B66*O_66(Jval),
'cubic_2':
    B40*O_40(Jval) + sympy.Rational(5,2)*B40*O_44(Jval) + \
    B60*O_60(Jval) + sympy.Rational(-21, 2)*B60*O_64(Jval)
}[lattice])

### Look into subspaces
def swap_matrix(dim:int, b1:int, b2:int):
    '''
    Create a matrix of rank `dim` that corresponds to swapping the basis vector
    number `b1` and `b2`.
    '''
    B = sympy.eye(dim)
    B[b1,b1] = B[b2,b2] = 0
    B[b1,b2] = B[b2,b1] = 1

    return B

def change_base(H: sympy.Matrix, subs: list):
    '''
    Determine the unitary operation that reorders the basis of the Hamiltonian
    according to the `subs` list.

    Returns: H, B
        H: sympy.Matrix
        Hamiltonian
        B: sympy.Matrix
        Matrix that changes the base
    '''
    dim = H.shape[0]
    B = sympy.eye(dim)
    for (b1, b2) in subs:
        B = B @ swap_matrix(dim, b1,b2)
    
    H = B.H @ H @ B

    return H, B
    
def to_blocks(H: sympy.Matrix, max_dim: int=99):
    '''
    Reorder the base to construct a block Hamiltonian.

    H
        Hamiltonian.
    max_dim
        Dimension up to which rearrange the Hamiltonian. Useful for debugging or diving into a problem.

    Returns: H, B
        H: sympy.Matrix
            Hamiltonian
        B: sympy.Matrix
            Matrix that changes the base
    '''
    dim = H.shape[0]
    B_cumultant = sympy.eye(dim)
    col, min_row = 0, 0
    for col in np.arange(min(dim, max_dim)):
        ss = np.where([not x.is_zero for x in H[min_row:, col]])[0]
        nominal_order = np.arange(len(ss))

        if not np.allclose(ss, nominal_order):
            subs = [(min_row+b1,min_row+b2) for b1,b2 in zip(nominal_order, ss)]
            H, B = change_base(H, subs)
            B_cumultant = B_cumultant @ B

        min_row += len(ss)

    return H, B_cumultant

H_untouched = H

sympy.pprint(H)
# H = change_base(H, [(1,2),(2,4),(3,6),(4,8),(5,10),(6,12)])
# H = change_base(H, [(1,4),(2,8),(4,7),(6,7)])
H, B = to_blocks(H)
# H = to_blocks(H)
# sympy.pprint(H)

basis = [freeionkets[b] for b in np.where(B.T)[1]]
starting_dim = 0 # This is not working anymore
for it,H_sub in enumerate(H.get_diag_blocks()):
    subspace_dim = H_sub.shape[0]

    if it==0 and False:
        print(f'Diagonalize Hamiltonian: {filenameprefix}')

        start_time = timer()
        P, D = H_sub.diagonalize()
        end_time = timer()

        diag_str = f'# Diagonalization time = {end_time-start_time:4f} s'
        print(diag_str)


        start_time = timer()
        Ds = sympy.nsimplify(D)
        Ps = sympy.nsimplify(P)
        end_time = timer()

        sim_str = f'# Simplification and saving time = {end_time-start_time:4f} s'
        print(sim_str)

        print(Ds)
        print(Ps)


    print(f'Subspace of dimension={subspace_dim}')
    print(f'Basis vectors:')
    print(basis[starting_dim:starting_dim+subspace_dim])
    print('Subspace Hamiltonian')
    sympy.pprint(H_sub)
    print()
    starting_dim += subspace_dim



sympy.pprint( B )

# Check if the change of basis didn't mess up the Hamiltonian
if not B@B.H==sympy.eye(B.shape[0]):
    raise ValueError(f'Change of base transformation is not hermitian')
if not H.is_hermitian:
    raise ValueError(f'Hamiltonian is not hermitian')



# Run calcualtions

if diagonalizeH:
    print(f'Diagonalize Hamiltonian: {filenameprefix}')

    start_time = timer()
    P, D = H.diagonalize()
    end_time = timer()

    diag_str = f'# Diagonalization time = {end_time-start_time:4f} s'
    print(diag_str)


    start_time = timer()
    Ds = sympy.nsimplify(D)
    Ps = sympy.nsimplify(P)
    end_time = timer()

    sim_str = f'# Simplification and saving time = {end_time-start_time:4f} s'
    print(sim_str)


### Latex
latex_begin = '''\\documentclass[8pt]{report}
\\usepackage[a4paper,margin=0.1in,landscape]{geometry}
\\usepackage{amsmath}
\\usepackage{graphicx}

\\begin{document}

'''
latex_end = '''
\\end{document}
'''

def eq_wrapper(text, resizebox=False):
    ret = f'\\begin{{math}}\n{text}\n\\end{{math}}'
    if resizebox:
        ret = '\\resizebox{0.98\\linewidth}{!}{%\n'+ret+'\n}'

    return ret

if exportLatex:
    hamiltonian_tex = sympy.latex(H)
    hamiltonian_un_tex = sympy.latex(H_untouched)

    # with open(f'{filenameprefix}.tex','w') as ff:
    with open(f'./Symbolic-output/test.tex','w') as ff:
        ff.write(latex_begin)
        
        ff.write('Hamiltonian:\n\n')
        ff.write(f'{eq_wrapper(hamiltonian_un_tex, resizebox=True)}\n\n')
        ff.write('Hamiltonian in blocks:\n\n')
        ff.write(f'{eq_wrapper(hamiltonian_tex, resizebox=True)}\n\n')

        basis = [freeionkets[b] for b in np.where(B.T)[1]]
        starting_dim = 0
        for H_sub in H.get_diag_blocks():
            subspace_dim = H_sub.shape[0]
            
            ff.write(f'Subspace of dimension={subspace_dim}\n\n')
            ff.write(f'Basis vectors: \n')
            str_basis = eq_wrapper(basis[starting_dim:starting_dim+subspace_dim]).replace("\'","")
            ff.write(f'{str_basis}\n\n')
            ff.write('Subspace Hamiltonian\n')
            ff.write(f'{eq_wrapper(sympy.latex(H_sub))}\n')
            ff.write(f'\n\n')
            starting_dim += subspace_dim
        
        if diagonalizeH:
            for it in range(int(2*Jval+1)):
                eigenvalue = eq_wrapper(sympy.latex(D[it,it]))
                eigenvector= eq_wrapper(sympy.latex(P[:,it]))
                ff.write(f'Eigenvalue and eigenvector {it+1}:\n\n')
                ff.write(f'{eigenvalue}\n\n')
                ff.write(f'{eigenvector}\n\n')
        
        ff.write(latex_end)


### Export results
if exportSympy and diagonalizeH:
    with open(f'{filenameprefix}-nsimplify.txt','w') as ff:
        ff.write(diag_str)
        ff.write('\n')
        ff.write(sim_str)
        ff.write('\n')
        ff.write('\n')

        
        for it in range(int(2*Jval+1)):
            ff.write(f'Eigenvalue and eigenvector {it+1}\n')
            ff.write(f'{Ds[it,it]}\n')
            ff.write(f'{Ps[:,it]}\n')
            ff.write(f'\n')
            
        ff.write(f'Diagonalized matrix D from H=U_dag D U\n')
        ff.write(f'{Ds}')
        ff.write(f'\n')
        ff.write(f'Rotation matrix U from H=U_dag D U\n')
        ff.write(f'{Ps}')