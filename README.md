# crysfipy-symbolic

Symbolic calculations in sympy aiming at determining ananlytical formulas for energy levels of the f-electron in different symmetry environments. Mathematical formalism is base on 

## To read
- https://en.wikipedia.org/wiki/Crystal_field_theory

## Sebastian
- Edit the `CEF-symbolic.py` script to define the J (`Jval`), symmetry (`Bij`) and the hamiltonian `H`.

## TASKS
1. [ ] Code in which `Bij` parameters go into Hamiltonian, depending on symmetry and implement into a switch case.
2. [ ] See if there is influence on the simplified formulas by changing floats in cefmatrices.py into sympy.Rationals and imaginary unit 'j' into sympy.I
3. [ ] Looking at https://ptable.com/#Electrons/OxidationStates determine J fo all lanthanides.

## Ideas
Ideas how to improve:
1. Calculate only eigenvalues -> thats actually super fast!!!
2. Start by reducing columns and rows -> Seem like sympy is smart enough to do that on its own

Ideas how to analyse
1. If the system is overdefined one can assume some parameters are small, and substitute them in equations to 0. This will produce a system with reduced number of parameters which will act as a initial set of parameters that can be later fit
2. Fin out about the physical boundaries on the parameters, like I think B_20>0
