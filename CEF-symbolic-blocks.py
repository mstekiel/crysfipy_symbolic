import sympy
import numpy as np

from cef_hamiltonian import CEF_Hamiltonian

from timeit import default_timer as timer


# Dont forget to edit these fields
Jval = 2.5
lattice='tetragonal_2'
# Ex.1 interesting set of 4D matrices
# Jval = 12
# lattice='hexagonal_2'
# Ex.2 5D matrix
# Jval = 12
# lattice='hexagonal_2'

diagonalizeH = False
investigateSubspaces = True

exportLatex = True
exportSympy = False

hline = '--------------------------------------------------------------------------------------'

###
# Calculations

H = CEF_Hamiltonian(lattice=lattice, Jval=Jval)
freeionkets = H.freeionkets


Jvalstr = f'{Jval:.1f}'.replace('.','p')
filenameprefix = f'./Symbolic-output/J_{Jvalstr}_{lattice}'



############################
# Diagonalization methods

def eig_values(H):
    dim = H.shape[0]
    # print(H[1])

    if dim==1:
        return H[0]
    elif dim==2:
        # [a, c]
        # [c, b]
        a, b, c = H[0], H[3], H[1]
        lam_1 = 1/2*(-sympy.sqrt(a**2 - 2*a*b + b**2 + 4*c**2) + a + b)
        lam_2 = 1/2*( sympy.sqrt(a**2 - 2*a*b + b**2 + 4*c**2) + a + b)
        return [lam_1.simplify(), lam_2.simplify()]
    elif dim==3:
        # [a, d, f]
        # [d, b, e]
        # [f, e, c]
        a, d, f = H[0], H[1], H[2]
        _, b, e = H[3], H[4], H[5]
        _, _, c = H[6], H[7], H[8]
        # P, D = H.diagonalize()
        # lam_1 = (2 a**3 - 3 a**2 b - 3 a**2 c + sympy.sqrt(4 (-a**2 + a b + a c - b**2 + b c - c**2 - 3 d**2 - 3 e**2 - 3 f**2)**3 + (2 a**3 - 3 a**2 b - 3 a**2 c - 3 a b**2 + 12 a b c - 3 a c**2 + 9 a d**2 - 18 a e**2 + 9 a f**2 + 2 b**3 - 3 b**2 c - 3 b c**2 + 9 b d**2 + 9 b e**2 - 18 b f**2 + 2 c**3 - 18 c d**2 + 9 c e**2 + 9 c f**2 + 54 d e f)**2) - 3 a b**2 + 12 a b c - 3 a c**2 + 9 a d**2 - 18 a e**2 + 9 a f**2 + 2 b**3 - 3 b**2 c - 3 b c**2 + 9 b d**2 + 9 b e**2 - 18 b f**2 + 2 c**3 - 18 c d**2 + 9 c e**2 + 9 c f**2 + 54 d e f)**(1/3)/(3 2**(1/3)) - (2**(1/3) (-a**2 + a b + a c - b**2 + b c - c**2 - 3 d**2 - 3 e**2 - 3 f**2))/(3 (2 a**3 - 3 a**2 b - 3 a**2 c + sqrt(4 (-a**2 + a b + a c - b**2 + b c - c**2 - 3 d**2 - 3 e**2 - 3 f**2)**3 + (2 a**3 - 3 a**2 b - 3 a**2 c - 3 a b**2 + 12 a b c - 3 a c**2 + 9 a d**2 - 18 a e**2 + 9 a f**2 + 2 b**3 - 3 b**2 c - 3 b c**2 + 9 b d**2 + 9 b e**2 - 18 b f**2 + 2 c**3 - 18 c d**2 + 9 c e**2 + 9 c f**2 + 54 d e f)**2) - 3 a b**2 + 12 a b c - 3 a c**2 + 9 a d**2 - 18 a e**2 + 9 a f**2 + 2 b**3 - 3 b**2 c - 3 b c**2 + 9 b d**2 + 9 b e**2 - 18 b f**2 + 2 c**3 - 18 c d**2 + 9 c e**2 + 9 c f**2 + 54 d e f)**(1/3)) + 1/3 (a + b + c)
        return [a, b, c]

############################
# Continue calculations

B = H.make_block_form()
sympy.pprint(H.original)
print('My own blocks')
sympy.pprint(H.in_blocks)
print(hline)

# print('Jordan form')
# P, J = H.original.jordan_form() # This looks interesting, but is slow for matrices with rank>=5
# sympy.pprint(P)
# sympy.pprint(J)


# H = to_blocks(H)
# sympy.pprint(H)

basis = [freeionkets[b] for b in np.where(B.T)[1]]
eigenvalues_all = {}
starting_dim = 0 # This is not working anymore
for it,H_sub in enumerate(H.in_blocks.get_diag_blocks()):
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

    print(hline)
    print(f'Subspace of dimension={subspace_dim}')
    print(f'Basis vectors:')
    print(basis[starting_dim:starting_dim+subspace_dim])
    print('Subspace Hamiltonian')
    sympy.pprint(H_sub)
    print('Subspace Hamiltonian eigenvalues')
    
    # for eig_value in eig_values(H_sub):
    #     sympy.pprint(eig_value)

    charpoly = H_sub.charpoly()
    # sympy.pprint(charpoly)
    # print(sympy.roots(charpoly))
    eigenvalues_subspace = sympy.roots(charpoly)

    for eigenvalue in eigenvalues_subspace:
        print(eigenvalue)
        eigenvectors = (H_sub - eigenvalue*sympy.eye(H_sub.rank())).nullspace()
        print(eigenvectors)

    # Fill the main dictionary
    for egv, degeneracy in eigenvalues_subspace.items():
        egv = egv.simplify()
        if egv in eigenvalues_all:
            eigenvalues_all[egv] += degeneracy
        else:
            eigenvalues_all[egv] = degeneracy

    starting_dim += subspace_dim

# Grouped eigenvalues
print(hline)
print(hline)
print('All eigenvalues (degeneration: value)')
for key, value in eigenvalues_all.items():
    print(value)
    print(sympy.factor(key) )
    # print(key.free_symbols)
    # B40 = list(key.free_symbols)[1]
    # B60 = list(key.free_symbols)[0]
    # B40 = sympy.symbols("B40")
    # print(key.subs(B60, 0) )

# Overall transformation matrix
# sympy.pprint( B )

# Check if the change of basis didn't mess up the Hamiltonian
if not B@B.H==sympy.eye(H.dimensionality):
    raise ValueError(f'Change of base transformation is not hermitian')
if not H.original.is_hermitian:
    raise ValueError(f'Hamiltonian (original) is not hermitian')
if not H.in_blocks.is_hermitian:
    raise ValueError(f'Hamiltonian (in block form) is not hermitian')



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
    hamiltonian_un_tex = sympy.latex(H.original)

    # with open(f'{filenameprefix}.tex','w') as ff:
    with open(f'./Symbolic-output/test.tex','w') as ff:
        ff.write(latex_begin)
        
        ff.write('Hamiltonian:\n\n')
        ff.write(f'{eq_wrapper(hamiltonian_un_tex, resizebox=True)}\n\n')
        ff.write('Hamiltonian in blocks:\n\n')
        ff.write(f'{eq_wrapper(hamiltonian_tex, resizebox=True)}\n\n')

        basis = [freeionkets[b] for b in np.where(B.T)[1]]
        starting_dim = 0
        for H_sub in H.original.get_diag_blocks():
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