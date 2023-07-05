import sympy
import numpy as np

from typing import Tuple

from  cef_matrices import *

class CEF_Hamiltonian():  
    def __init__(self, lattice: str, Jval: float):
        self.dimensionality = int(2*Jval+1)

        # Symbols for Hamiltonian
        B20, B22, B2m2 = sympy.symbols('B20, B22, B2m2', real=True)
        B40, B42, B4m2, B43, B4m3, B44, B4m4 = sympy.symbols('B40, B42, B4m2, B43, B4m3, B44, B4m4', real=True)
        B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6 = sympy.symbols('B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6', real=True)

        # Symmetry allowed Bnm according to McPhase manual
        # https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html
        Hamiltonian = sympy.Matrix({
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
            B40*O_40(Jval) + sympy.Integer(5)*B40*O_44(Jval) + \
            B60*O_60(Jval) - B66*O_62(Jval) + sympy.Integer(-21)*B60*O_64(Jval) + B66*O_66(Jval),
        'cubic_2':
            B40*O_40(Jval) + sympy.Integer(5)*B40*O_44(Jval) + \
            B60*O_60(Jval) + sympy.Integer(-21)*B60*O_64(Jval)
        }[lattice])

        self.original = Hamiltonian
        self.in_blocks = Hamiltonian

        # Prepare a list of kets that form the basis
        if Jval%1==0:
            freeionkets = [f'|{int(x)}>' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]
        elif Jval%1==0.5:
            freeionkets = [f'|{int(2*x):d}/2>' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]

        self.freeionkets = freeionkets


    # Methodes for making the Jordan form of the Hamiltonian
        
    ### Look into subspaces
    def _swap_matrix(self, dim:int, b1:int, b2:int) -> sympy.Matrix:
        '''
        Create a matrix of rank `dim` that corresponds to swapping the basis vector
        number `b1` and `b2`.
        '''
        B = sympy.eye(dim)
        B[b1,b1] = B[b2,b2] = 0
        B[b1,b2] = B[b2,b1] = 1

        return B

    def _change_base_matrix(self, subs: list) -> Tuple[sympy.Matrix, sympy.Matrix]:
        '''
        Determine the unitary operation that reorders the basis of the Hamiltonian
        according to the `subs` list.
        
        subs: List[Tuple[int,int]]
            List of substitutions between base vectors.

        Returns: B
            B: sympy.Matrix
            Matrix that changes the base
        '''
        dim = self.dimensionality
        B = sympy.eye(dim)
        for (b1, b2) in subs:
            B = B @ self._swap_matrix(dim, b1,b2)
        
        return B
        
    def make_block_form(self, max_dim: int=99) -> Tuple[sympy.Matrix, sympy.Matrix]:
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
        dim = self.dimensionality
        B_cumultant = sympy.eye(dim)
        col, min_row = 0, 0
        for col in np.arange(min(dim, max_dim)):
            # Find entries in a column that are non zero, and output only row numbers
            block_order = np.where([not x.is_zero for x in self.in_blocks[min_row:, col]])[0]
            # Make a list of nominal order to aid constructing substitutions
            nominal_order = np.arange(len(block_order))

            # Sometimes the block order is the same as nominal, and its better not to change it
            if not np.allclose(block_order, nominal_order):
                subs = [(min_row+b1,min_row+b2) for b1,b2 in zip(nominal_order, block_order)]
                B = self._change_base_matrix(subs)
                self.in_blocks = B.H @ self.in_blocks @ B
                B_cumultant = B_cumultant @ B

            min_row += len(block_order)

        return B_cumultant
