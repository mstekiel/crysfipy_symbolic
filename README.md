# crysfipy-symbolic

Symbolic calculations in sympy aiming at determining ananlytical formulas for energy levels of the f-electron in different symmetry environments.

## To read
- https://en.wikipedia.org/wiki/Crystal_field_theory

## TODO
- [ ] Take a look at the Hamiltonian for all symmetries and J up to 5 and determine what kind of blocks are present.
- [ ] Based on how the blocks look like, implement the solutions for eigenvalues for specific matrices.
- [X] Code in which `Bij` parameters go into Hamiltonian, depending on symmetry.
  - [X] please double-check my new implementation with www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html webpage (corrected, see commit history)
- [X] Run calculations for various symmetries and check what is the highes J for which the calculations runs fast (let's say up to 10 seconds).
- [X] See if there is influence on the simplified formulas by changing floats in cefmatrices.py into sympy.Rationals.

## Ideas
Ideas how to improve:
- [X] Calculate only eigenvalues -> thats actually super fast!!!
- [X] Start by reducing columns and rows -> Seem like sympy is smart enough to do that on its own

Ideas how to analyse
- If the system is overdefined one can assume some parameters are small, and substitute them in equations to 0. This will produce a system with reduced number of parameters which will act as a initial set of parameters that can be later fit
- Find out about the physical boundaries on the parameters, like I think B_20>0
