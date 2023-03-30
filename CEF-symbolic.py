import sympy
import numpy as np

from timeit import default_timer as timer

from cefmatrices import *

sympy.init_printing(use_unicode=False)
B20, B22, B2m2 = sympy.symbols('B20, B22, B2m2', real=True)
B40, B42, B4m2, B43, B4m3, B44, B4m4 = sympy.symbols('B40, B42, B4m2, B43, B4m3, B44, B4m4', real=True)
B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6 = sympy.symbols('B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6', real=True)


# Dont forget to edit these fields
Jval = 2
lattice='monoclinic'


###
# Calculations

Jvalstr = f'{Jval:.1f}'.replace('.','p')
filenameprefix = f'./Symbolic-output/J_{Jvalstr}_{lattice}'

# Symmetry allowed Bnm according to McPhase manual
# https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html
H = sympy.Matrix({
'monoclinic':
    B20*O_20(Jval) + B22*O_22(Jval) + B2m2*O_2m2(Jval) + \
    B40*O_40(Jval) + B42*O_42(Jval) + B4m2*O_4m2(Jval) + B44*O_44(Jval) + B4m4*O_4m4(Jval) + \
    B60*O_60(Jval) + B62*O_62(Jval) + B6m2*O_6m2(Jval) + \
    B64*O_64(Jval) + B6m4*O_6m4(Jval) + B66*O_66(Jval) + B6m6*O_6m6(Jval),
'rhombic':
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
    B40*O_40(Jval) + sympy.Rational(5,2)*B40*O_40(Jval) + \
    B60*O_60(Jval) + -1*B66*O_66(Jval) + sympy.Rational(-21, 2)*B60*O_60(Jval) + B66*O_66(Jval),
'cubic_2':
    B40*O_40(Jval) + sympy.Rational(5,2)*B40*O_40(Jval) + \
    B60*O_60(Jval) + sympy.Rational(-21, 2)*B60*O_60(Jval)
}[lattice])

sympy.pprint( H )
H_dag = sympy.transpose(sympy.conjugate(H))
if not H==H_dag:
    sympy.pprint(H_dag)
    raise ValueError(f'Hamiltonian is not hermitian')

diagonalize = True
exportLatex = False
exportSympy = True
exportMD = True



# Run calcualtions

print(f'Run calculations for Hamiltonian: {filenameprefix}')

hamiltonian = sympy.latex(H)

start_time = timer()

if diagonalize:
    P, D = H.diagonalize()

    eigenvalues = D
    eigenvectors = P

    # print(D)
    # print()
    # print(P)
else:
    eigenvalues = '?'
    eigenvectors = '?'


    
end_time = timer()
diag_str = f'Diagonalization time = {end_time-start_time} s'
print(diag_str)


latex_begin = '''//documentclass[8pt]{report}
//usepackage[a4paper,margin=0.1in,landscape]{geometry}
//usepackage{amsmath}

//begin{document}

'''
latex_end = '''
//end{document}
'''

def eq_wrapper(text):
    return f'//begin{{math}}\n{text}\n//end{{math}}'


start_time = timer()

if exportLatex:
    with open(f'{filenameprefix}.tex','w') as ff:
        ff.write(latex_begin)
        
        ff.write('Hamiltonian:\n\n')
        ff.write(f'{eq_wrapper(hamiltonian)}\n\n')
        
        for it in range(int(2*Jval+1)):
            eigenvalue = eq_wrapper(sympy.latex(D[it,it]))
            eigenvector= eq_wrapper(sympy.latex(P[:,it]))
            ff.write(f'Eigenvalue and eigenvector {it+1}:\n\n')
            ff.write(f'{eigenvalue}\n\n')
            ff.write(f'{eigenvector}\n\n')
        
        ff.write(latex_end)
        
if exportSympy:
    with open(f'{filenameprefix}-nsimplify.txt','w') as ff:
        Ds = sympy.nsimplify(D)
        Ps = sympy.nsimplify(P)
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
            
        
end_time = timer()
sim_str = f'Simplification and saving time = {end_time-start_time} s'
print(sim_str)

if exportMD:
    with open(f'{filenameprefix}.md','w') as ff:
        ff.write(diag_str)
        ff.write('\n')
        ff.write(sim_str)
