# crysfipy_-_symbolic

Symbolic calculations in sympy aiming at determining ananlytical formulas for energy levels of the f-electron in different symmetry environments.

## To read
- https://en.wikipedia.org/wiki/Crystal_field_theory
- http://cbpfindex.cbpf.br/publication_pdfs/NF01594.2011_06_07_15_09_34.pdf Brazilians did obtain analytical results for the cubic case with limited number of J values. They are different from ours. We need to check numerically if our definitions are ok.

## To check
- https://www.sciencedirect.com/science/article/pii/0022369762901920 apparently a seminal work work the 'cubic_2' case.

## TODO
- [ ] Simplify the formulas, by looking at characterstic form, e.g. (a+b)/(a+c) `find repeating substrings`, wrapping imaginary unit into sqrt etc. -> look into the `sympy.cse` function.
- [ ] Determine coefficients that occur only in squared form, as the sign of those will be impossible to determine from dE. However, there are some cases that also with linear terms the sign can be ambigous due to quantum phase.
- [ ] Determine conditions on parameters that need to be fulfilled so that energies are real -> this should be automatic... -> check that eigenvalues are real
- [X] Implement yaml dumper -> using automatic protocol from `pyaml` library.
- [X] Test time needed to run calculations -> results in `test` folder.
- [X] Proper architecture of the dependencies. 
- [X] Mapping between lattice-Schoenflies-HermannMauguin noteation for point group symmetry.
- [X] Construct the characteristic polynomial and investigate whether it can be divided by other polynomial. -> This works marvelously !!!
  - [X] Take a look at the eigenvectors
- [X] Take a look at the Hamiltonian for all symmetries and J up to 5 and determine what kind of blocks are present.
- [X] Based on how the blocks look like, implement the solutions for eigenvalues for specific matrices -> depracted, getting roots of characteristic polynomial works well.
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


## Architecture
Hamiltonian
- has base, and matrix as defining properties
- diagonalization 
	- costly process, should be done only once by function call
 	- calculating eigenvalues is way faster than calculating eigenvectors, so separate function should be used to call each procedure
	- eigenvalues can be degenerated


CEF Hamiltonian
- is a Hamiltonian with additional properties:
	- has `symmetry` and `J` that define the base and matrix of Hamiltonian
	- can always? be reduced to block form
- block form
	- each block is a Hamiltonian
	- eigenvalues from different blocks may be the same, in fact, there's always repeating eigenvaluesin case of J = half-integer (Kramers theorem)
	- how to store those repeated eigenvalues?
- diagonalization should be done only in the block form, as it's too costly to diagonalize full matrix
- how to store the eigenvalues/eigenvectors? override the inherited fields or make new ones?
	
	
### Icing
- nice printing
- latex printing
- yaml dimping
- numerical evaluations