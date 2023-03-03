# crysfipy-symbolic

Symbolic calculations in sympy aiming at determining ananlytical formulas for energy levels of the f-electron in different symmetry environments. Mathematical formalism is base on 

## To read
- https://en.wikipedia.org/wiki/Crystal_field_theory

## TODO
1. See if there is influence on the simplified formulas by changing floats in cefmatrices.py into sympy.Rationals.

## Ideas
Ideas how to improve:
1. Calculate only eigenvalues -> thats actually super fast!!!
2. Start by reducing columns and rows -> Seem like sympy is smart enough to do that on its own

Ideas how to analyse
1. If the system is overdefined one can assume some parameters are small, and substitute them in equations to 0. This will produce a system with reduced number of parameters which will act as a initial set of parameters that can be later fit
2. Fin out about the physical boundaries on the parameters, like I think B_20>0