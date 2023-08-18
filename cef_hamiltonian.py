import sympy
import numpy as np
from copy import deepcopy

from typing import List, Dict
from sympy.matrices import Matrix, ShapeError
from cef_matrices import *

from sympy import pretty

class Hamiltonian():
    '''
    Class representing a Hamiltonian.

    Fields
    ------
    base: List[str]
        Base of the Hamiltonian.
    matrix: sympy.Matrix
        Matrix representation of the hamiltonian in `self.base`
    '''

    def __init__(self, base: List[str], matrix: Matrix):
        if not matrix.shape[0] == matrix.shape[1]:
            raise ShapeError(f'Hamiltonian matrix is not square. Defined dimensions={matrix.shape}')
        
        if not len(base) == matrix.shape[1]:
            raise ShapeError(f'Base ({len(base)},) does not align with matrix representation {matrix.shape}')

        self.base = base
        self.matrix = matrix

        self.eigenvalues = None
        self.eigenvectors = None

    def __str__(self):
        ret = 'Hamiltonian base\n'
        ret += str(self.base) + '\n'
        ret += 'Hamiltonian matrix representation\n'
        ret += pretty(self.matrix)

        return ret
    
    def cal_eigenvalues(self):
        '''
        Determine the eigenvalues of the Hamiltonian by constructing the
        characteristic polynomial and finding its roots.
        '''
        self.eigenvalues = sympy.roots(self.matrix.charpoly())
        return self.eigenvalues

    def sort_base(self):
        '''
        Sort the base of the Hamiltonian according to the numerical representation
        of the base.

        Parameters
        ----------        
        new_order: List[int]
            List of unique integers between 0 and dimensionality of Hamiltonian.
            Reorders the base

        Update fields
        -------------
        base
        matrix
        '''

        # Sort descending order of the number representation
        sorted_order = np.argsort([eval(x) for x in self.base])[::-1]
        B = self.reorder_base(sorted_order)


    def reorder_base(self, new_order: List[int]) -> Matrix:
        '''
        Determine the unitary operation that reorders the basis of the Hamiltonian
        according to the `new_order` list.

        Parameters
        ----------        
        new_order: List[int]
            List of unique integers between 0 and dimensionality of Hamiltonian
            that reorders the base.

        Returns
        -------
        B: sympy.Matrix
            Matrix that reorders the base
        '''
        B = sympy.zeros(len(new_order))
        for b1, b2 in enumerate(new_order):
            B[b1,b2] = 1

        self.matrix = B @ self.matrix @ B.T
        self.base = [self.base[i] for i in new_order]

        if self.eigenvalues:
            self.eigenvalues = [self.eigenvalues[i] for i in new_order]

        if self.eigenvectors:
            self.eigenvectors = [self.eigenvectors[i] for i in new_order]

        
        return B


class CEF_Hamiltonian(Hamiltonian):
    '''
    Class representing the Hamiltonian of an ion in crystal field.

    Fields
    ------
    base: List[str]
        List of ordered kets forming the basis.
    matrix: sympy.Matrix
        Hamiltonian matrix in the original form i.e., basis given by `freeionkets`.
    
    subs: list[Hamiltonian]
        Subhamiltonians extracted from `self.make_block_form`.
    '''
    def __init__(self, symmetry: str, Jval: float):
        # Prepare a list of kets that form the basis
        if Jval%1==0:
            freeionkets = [f'{int(x)}' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]
        elif Jval%1==0.5:
            freeionkets = [f'{int(2*x):d}/2' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]

        # Symbols for Hamiltonian
        B20, B22, B2m2 = sympy.symbols('B20, B22, B2m2', real=True)
        B40, B42, B4m2, B43, B4m3, B44, B4m4 = sympy.symbols('B40, B42, B4m2, B43, B4m3, B44, B4m4', real=True)
        B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6 = sympy.symbols('B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6', real=True)

        # Symmetry allowed Bnm according to McPhase manual
        # https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html
        matrix_representation = sympy.Matrix({
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
        }[symmetry])

        # Constructor of the parent class
        super().__init__(base=freeionkets, matrix=matrix_representation)

        # Other class specific fields
        self.subs = None




    # Methods for making the Jordan form of the Hamiltonian
    def make_block_form(self) -> List[Hamiltonian]:
        '''
        Reorder the base to construct a block Hamiltonian.

        Updated fields
        --------------
        subs

        Returns
        -------
        blocks
            Subhamiltonians obtained from Jordan form the block form of the Hamiltonian
        '''

        from sympy import pprint

        blocks = []

        H_red = deepcopy(self)
        # H_red = deepcopy(self.original)
        col = 0
        while col < len(self.base):
            # Find entries in a row that are non zero and reorder
            non_zero_entries = list(np.where([not x.is_zero for x in H_red.matrix[0, :]])[0])
            zero_entries = list(np.where([x.is_zero for x in H_red.matrix[0, :]])[0])
            new_order = non_zero_entries + zero_entries
            sub_rank = len(non_zero_entries)

            H_red.reorder_base(new_order)

            # Extract sub Hamiltonian
            H_sub = Hamiltonian(
                base = H_red.base[:sub_rank],
                matrix = H_red.matrix[:sub_rank, :sub_rank]
            )
            H_sub.sort_base()
            blocks.append(H_sub)

            # Reduce the Hamiltonian for the next iteration
            H_red = Hamiltonian(
                base = H_red.base[sub_rank:],
                matrix = H_red.matrix[sub_rank:, sub_rank:]
            )

            col += sub_rank

        self.subs = blocks

        return blocks