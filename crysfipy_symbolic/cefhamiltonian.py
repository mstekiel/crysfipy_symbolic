# Tools
import numpy as np
from copy import deepcopy
import warnings
import yaml

# Typesetting
from typing import List, Dict, Tuple, Any
from sympy import Expr as expr

# Sympy
import sympy
from sympy import pretty
from sympy.matrices import Matrix, ShapeError

# Internal
from .cefmatrices import *
from .hamiltonian import Hamiltonian
from .data_containers import sympy_expression_mapping, point_group_synonyms


class CEF_Hamiltonian(Hamiltonian):
    '''
    Class representing the Hamiltonian of an ion in crystal field.

    Fields
    ------
    Jval: float
        J quantum number.
    base: List[str]
        List of ordered kets forming the basis.
    matrix: sympy.Matrix
        Hamiltonian matrix in the original form i.e., basis given by `freeionkets`.
    symmetry: Tuple[str, str, str]
        Point group names: (Hermann-Mauguin, Schoenflies, lattice).
    
    subs: list[Hamiltonian]
        Subhamiltonians extracted from `self.make_block_form`.
    eigenvalues: list[tuple[expr, int, int]]
        List of eigenvalues with formula, degeneracy and subspace origin index.
    eigenvectors: list[tuple[expr, int, int]]
        List of eigenvectors with formula, index of eigenvalue from which it originates, and subspace origin index.
    '''
    def __init__(self, symmetry: str, Jval: float):
        '''
        Construct crystal field Hamiltonian with `Jval` quantum number in an environment with specific `symmetry`.
        '''

        if not np.isclose( 2*Jval, int(2*Jval)):
            raise ValueError(f'J={Jval} is not half-integer.')

        self.Jval = float(Jval)

        # Prepare a list of kets that form the basis
        if int(2*self.Jval) % 2:
            freeionkets = [f'{int(2*x):d}/2' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]
        else:
            freeionkets = [f'{int(x)}' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]


        # Symmetry allowed Bnm according to McPhase manual
        # https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html
        matrix_representation = self.H_from_point_group(symmetry, Jval)
        self.symmetry =  point_group_synonyms(symmetry)[0] # take the first suggested synonym

        # Constructor of the parent class
        super().__init__(base=freeionkets, matrix=matrix_representation)

        # Other class specific fields
        self.subs = []

    def H_from_point_group(self, point_group_name: str, Jval: float) -> Matrix:
        '''
        Construct Hamiltonian matrix based on the point group symbol.

        Symmetry allowed Bnm according to McPhase manual
        https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html

        Parameters:
        -----------
        point_group_name
            String representing the point group in lattice, Schoenflies, Hermann-Maguire notation.
            In case of ambiguity, the first matching point group symmetry is taken.

        Returns:
        --------
            Matrix representation of the Hamiltonian
        '''

        hm, scheonflies, lattice, id = point_group_synonyms(point_group_name)[0]
        symmetry = (hm, scheonflies, lattice)

        # Symbols for Hamiltonian
        B20, B22, B2m2 = sympy.symbols('B20, B22, B2m2', real=True)
        B40, B42, B4m2, B43, B4m3, B44, B4m4 = sympy.symbols('B40, B42, B4m2, B43, B4m3, B44, B4m4', real=True)
        B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6 = sympy.symbols('B60, B62, B6m2, B63, B6m3, B64, B6m4, B66, B6m6', real=True)

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
        }[lattice])

        return matrix_representation

    ##################################################################################
    # Diagonalization
    def determine_eigenvalues(self, from_blocks: bool=True) -> List[Tuple[expr, int, int]]:
        '''
        Determine eigenvalues of the Hamlitonian
        '''
        if not from_blocks:
            return super().determine_eigenvalues()
        
        if not self.subs:
            self.make_block_form()

        # This dict will store eigenvalues, their degeneracy and subspace origin
        # eigenvalue_string: [subID1, ...]
        eigenvalues_dict = {}
        for n,Hsub in enumerate(self.subs):
            for eigenvalue, degeneracy in Hsub.determine_eigenvalues():
                if eigenvalue in eigenvalues_dict:
                    eigenvalues_dict[eigenvalue].extend([n]*degeneracy)
                else:
                    eigenvalues_dict[eigenvalue] = [n]*degeneracy

        self.eigenvalues = [(ev, len(sub_ids), sub_ids) for ev,sub_ids in eigenvalues_dict.items()]

        return self.eigenvalues
    
    def determine_eigenvectors(self, from_blocks: bool=True) -> List[Tuple[expr, int, int]]:
        '''
        Updated fields
        --------------

        eigenvectors
        '''
        if not from_blocks:
            return super().determine_eigenvectors()
        
        if not self.subs:
            self.make_block_form()

        if not self.eigenvalues:
            self.determine_eigenvalues()

        self.eigenvectors = []
        for subspace_origin, Hsub in enumerate(self.subs):
            subspace_base = Hsub.base
            for subspace_eigenvector, eigenvalue_origin_in_subspace in Hsub.determine_eigenvectors():
                # Extend the base
                eigenvector = sympy.matrices.zeros(1, len(self.base))
                for vi, bi in zip(subspace_eigenvector, subspace_base):
                    eigenvector[self.base.index(bi)] = vi

                # Locate the index of the eigenvalue in whole Hamiltonian
                eigenvalue = Hsub.eigenvalues[eigenvalue_origin_in_subspace][0]
                eigenvalue_origin = [eval[0] for eval in self.eigenvalues].index(eigenvalue)

                self.eigenvectors.append( (eigenvector, eigenvalue_origin, subspace_origin) )

        return self.eigenvectors


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

        # Loop through all rows, find non-zero elements in each,
        # see what subspaces they make by grouping indices which
        # construct the subspace
        blocks_ids = []
        for row in range(self.matrix.shape[0]):
            id_non_zero_entries = set([n for n,x in enumerate(self.matrix[row,:]) if not x.is_zero])

            # Find block which makes a subspace with current row entries
            make_new_block = True
            for block_ids in blocks_ids:
                if id_non_zero_entries.intersection(block_ids):
                    block_ids.update(id_non_zero_entries)
                    make_new_block = False
                    break

            # If the new entries don't fit in previous blocks
            # make a new one
            if make_new_block:
                blocks_ids.append(id_non_zero_entries)

        # Extract blocks
        blocks = []
        for block_ids in blocks_ids:
            it = list(block_ids)
            H_sub = Hamiltonian(
                base = [self.base[n] for n in it],
                matrix = self.matrix[it, :][:, it]
            )
            H_sub._sort_base()
            blocks.append(H_sub)

        # Save and return the results
        self.subs = blocks
        return blocks
    
    ############################################################################
    # Printing

    def to_dict(self) -> Dict[str, Any]:
        '''
        Represent as dictionary
        '''

        H_dict = {}

        H_dict['base'] = self.base
        H_dict['J'] = self.Jval
        H_dict['symmetry'] = str(self.symmetry)
        H_dict['free_pars'] = [str(x) for x in self.matrix.free_symbols]

        # Matrix -> 2x2 list
        H_dict['matrix'] = []
        dim = self.matrix.rank()
        for n in range(dim):
            row = self.matrix[n*dim: (n+1)*dim]
            H_dict['matrix'].append([str(aij) for aij in row])

        # Eigenvalues -> (eigenvalue, degeneracy subspace_origins)
        H_dict['eigenvalues'] = {}
        for n, (eval, deg, subspace_origins) in enumerate(self.eigenvalues):
            H_dict['eigenvalues'][n] = dict(formula=str(eval), degeneracy=deg, subspace_origins=subspace_origins)

        # Eigenvectors ->  (evec, eigenvalue_origin, subspace_origin) 
        H_dict['eigenvectors'] = {}
        for n, (evec, eigenvalue_origin, subspace_origin) in enumerate(self.eigenvectors):
            eigenvector = [str(u) for u in evec]
            H_dict['eigenvectors'][n] = dict(formula=eigenvector, subspace_origin=subspace_origin, eigenvalue_origin=eigenvalue_origin)

        H_dict['subspaces'] = {}
        for n, Hsub in enumerate(self.subs):
            H_dict['subspaces'][n] = Hsub.to_dict()

        return H_dict


    def to_latex_doc(self) -> str:
        '''
        Parse the fields that need to be printed into latex.
        '''
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

        _eq_wrapper = lambda text: f'\\begin{{math}}\n{text}\n\\end{{math}}'
        latex = sympy.printing.latex

        ret_latex = latex_begin
        ret_latex += _eq_wrapper(latex(self.matrix))

        for n_sub, Hsub in enumerate(self.subs):
            ret_latex += f'\n\n'
            ret_latex += f'Subspace {n_sub+1}\n\n'
            ret_latex += 'Base: ' + ', '.join([f'$\\vert {x} \\rangle$' for x in Hsub.base]) + '\n\n'
            ret_latex += _eq_wrapper(f'H_{{{n_sub+1}}} = \n' + latex(Hsub.matrix))
            ret_latex += f'\n\n'


            Hsub.determine_eigenvalues()
            for n_ev,ev in enumerate(Hsub.eigenvalues):
                ev_expr, ev_deg = ev
                ret_latex += _eq_wrapper(f'\\lambda_{{{n_sub+1},{n_ev+1}}}^{{({ev_deg})}} = \n' + latex(ev_expr))
                ret_latex += f'\n\n'


        ret_latex += latex_end

        return ret_latex
    
################################################################
# For quick testing
if __name__ == '__main__':
    Jval = 2.5
    symmetry = '4mm'

    H = CEF_Hamiltonian(symmetry=symmetry, Jval=Jval)
    H.determine_eigenvectors()
    print(H.eigenvectors)
    sympy.matrices.zeros(*H.eigenvectors[0][0].shape)