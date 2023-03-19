import sympy
import numpy as np

from timeit import default_timer as timer

from cefmatrices import *

sympy.init_printing(use_unicode=False)
B20, B22, B40, B42, B43, B44, B60, B62, B63, B64, B66 = sympy.symbols('B20, B22, B40, B42, B43, B44, B60, B62, B63, B64, B66')


# Dont forget to edit these fields
Jval = sympy.Rational(3,2)
lattice='hexagonal_1' #can vary

Jvalstr = f'{Jval:.1f}'.replace('.','p')
filenameprefix = f'./Symbolic-output/J_{Jvalstr}_lattice'

match lattice:
	case 'triclinic':
		H = sympy.Matrix(B20 * O_20(Jval) + B22 * O_22(Jval) + B40 * O_40(Jval) + B42 * O_42(Jval) + B43 * O_43(Jval) + B44 * O_44(Jval) + B60 * O_60(Jval) + B62 * O_62(Jval) + B63 * O_63(Jval) + B64 * O_64(Jval) + B66 * O_66(Jval))
	case 'monoclinic':
		H = sympy.Matrix(B20 * O_20(Jval) + B22 * O_22(Jval) + B40 * O_40(Jval) + B42 * O_42(Jval) + B44 * O_44(Jval) + B60 * O_60(Jval) + B62 * O_62(Jval) + B64 * O_64(Jval))
	case 'rhombic':
		H = sympy.Matrix(B20 * O_20(Jval) + B22 * O_22(Jval) + B40 * O_40(Jval) + B42 * O_42(Jval) + B44 * O_44(Jval) + B60 * O_60(Jval) + B62 * O_62(Jval) + B64 * O_64(Jval))
	case 'tetragonal_1':
		H = sympy.Matrix(B20 * O_20(Jval) + B40 * O_40(Jval) + B44 * O_44(Jval) + B60 * O_60(Jval) + B64 * O_64(Jval))
	case 'tetragonal_2':
		H = sympy.Matrix(B20 * O_20(Jval) + B40 * O_40(Jval) + B44 * O_44(Jval) + B60 * O_60(Jval) + B64 * O_64(Jval))
	case 'trigonal_1':
		H = sympy.Matrix(B20 * O_20(Jval) + B40 * O_40(Jval) + B43 * O_43(Jval) + B60 * O_60(Jval) + B63 * O_63(Jval) + B66 * O_66(Jval))
	case 'trigonal_2':
		H = sympy.Matrix(B20 * O_20(Jval) + B40 * O_40(Jval) + B43 * O_43(Jval) + B60 * O_60(Jval) + B63 * O_63(Jval) + B66 * O_66(Jval))
	case 'hexagonal_1':
		H = sympy.Matrix(B20 * O_20(Jval) + B40 * O_40(Jval) + B60 * O_60(Jval) + B66 * O_66(Jval))
	case 'hexagonal_2':
		H = sympy.Matrix(B20 * O_20(Jval) + B40 * O_40(Jval) + B60 * O_60(Jval) + B66 * O_66(Jval))
	case 'cubic_1':
		H = sympy.Matrix(B40 * O_40(Jval) + sympy.Rational(5,2) * B40 * O_40(Jval) + B60 * O_60(Jval) + -1 * B66 * O_66(Jval) + sympy.Rational(-21, 2) * B60 * O_60(Jval) + B66 * O_66(Jval))
	case 'cubic_2':
		H = sympy.Matrix(B40 * O_40(Jval) + sympy.Rational(5,2) * B40 * O_40(Jval) + B60 * O_60(Jval) + sympy.Rational(-21, 2) * B60 * O_60(Jval))

print( H )


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
            
        ff.write(f'{Ds}')
        ff.write(f'\n')
        ff.write(f'{Ps}')
            
        
end_time = timer()
sim_str = f'Simplification and saving time = {end_time-start_time} s'
print(sim_str)

if exportMD:
    with open(f'{filenameprefix}.md','w') as ff:
        ff.write(diag_str)
        ff.write('\n')
        ff.write(sim_str)
