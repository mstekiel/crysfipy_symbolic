import sympy
import sympy.physics.quantum.spin as ssympy
from sympy.physics.quantum import Operator, represent, Ket, qapply

import numpy as np

from timeit import default_timer as timer

from cefmatrices import *

sympy.init_printing(use_unicode=False)
B20, B22, B40, B42, B44, B60, B62, B64, B66 = sympy.symbols('B20, B22, B40, B42, B44, B60, B62, B64, B66')
m = sympy.symbols('m')


# Symbolic analysis of the crystal field effects in CeAuAl3

# Dont forget to edit these fields
Jval = sympy.Rational(5,2)
Jvalstr = f'{Jval:.1f}'.replace('.','p')
filenameprefix = f'C:/Users/Stekiel/Documents/GitHub/mikibox/mikibox/crysfipy/symbolic/Symbolic-output/J_{Jvalstr}_B20_B22'

J2 = ssympy.J2Op()
Jplus = ssympy.JplusOp()
Jminus = ssympy.JminusOp()
Jz = ssympy.JzOp()

Jz_state = ssympy.JzKet(Jval, m)
Jz_other = m*ssympy.JzKet(Jval, m)*ssympy.JzBra(Jval, m)

print(Jz)
print(qapply(Jz*ssympy.JzKet(Jval, m)))

# O_20 = 3*Jz*Jz - J2
# O_22 = (Jplus**2 + Jminus**2)/2

# H = B20 * O_20 + B22 * O_22

# print(H)
# print(represent(H, basis=ssympy.JzKet))


#H = sympy.Matrix(   B20 * O_20(Jval) + B22 * O_22(Jval)) # + B40 * O_40(Jval)) + B42 * O_42(Jval) + B44 * O_44(Jval))
