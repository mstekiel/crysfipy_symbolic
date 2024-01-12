import numpy as np
from copy import deepcopy
import warnings
import yaml

from typing import List, Dict, Tuple

import sympy
from sympy import pretty
from sympy.matrices import Matrix, ShapeError

from .cefmatrices import *
from .data_containers import sympy_expression_mapping, point_group_synonyms


class Hamiltonian():
    '''
    Class representing a Hamiltonian.

    Fields
    ------
    base: List[str]
        Base of the Hamiltonian.
    matrix: sympy.Matrix
        Matrix representation of the hamiltonian in `self.base`
    eigenvalues: List[Tuple[expr, int]]
        Lsit of eigenvalues, where each entry contains
        eigenvalue's formula and degeneracy.
    eigenvectors: List[List[expr]]
        List of eigenvectors in order of `eigenvalues`.

    '''

    def __init__(self, base: List[str], matrix: Matrix):
        '''
        Class constructor.

        Requires `base` and `matrix` fields to fill, sets others to `None`.
        Checks the dimensionality of the Hamiltonian by ensuring equal dimensions 
        of `base` and `matrix`.
        '''
        if not matrix.shape[0] == matrix.shape[1]:
            raise ShapeError(f'Hamiltonian matrix is not square. Defined dimensions={matrix.shape}')
        
        if not len(base) == matrix.shape[1]:
            raise ShapeError(f'Base ({len(base)},) does not align with matrix representation {matrix.shape}')

        self.base = base
        self.matrix = matrix

        self.eigenvalues = []
        self.eigenvectors = []

        self._translate_table = sympy_expression_mapping()

    ##############################################################################33
    # Printing and representing
    def __str__(self):
        ret = 'Hamiltonian base\n'
        ret += str(self.base) + '\n'
        ret += 'Hamiltonian matrix representation\n'
        ret += pretty(self.matrix)

        return ret
    
    def to_dict(self):
        '''
        Represent as dictionary
        '''

        H_dict = {}

        # Base
        H_dict['base'] = self.base

        # Free parameters
        H_dict['free_pars'] = [str(x) for x in self.matrix.free_symbols]

        # Matrix
        H_dict['matrix'] = []
        dim = self.matrix.rank()
        for n in range(dim):
            row = self.matrix[n*dim: (n+1)*dim]
            H_dict['matrix'].append([str(aij) for aij in row])

        # Eigenvalues
        H_dict['eigenvalues'] = {}
        for n, eval in enumerate(self.eigenvalues):
            H_dict['eigenvalues'][n] = dict(formula=str(eval[0]), degeneracy=eval[1])

        # Eigenvectors
        H_dict['eigenvectors'] = {}
        for n, (evec, eval_origin) in enumerate(self.eigenvectors):
            eigenvector = [str(u) for u in evec]
            H_dict['eigenvectors'][n] = dict(formula=eigenvector, eigenvalue_origin=eval_origin)

        return H_dict
    
    def to_yaml(self):
        '''
        Convert the Hamiltonian to YAML format.

        Works by calling the `self.to_dict()` function.
        '''
        return yaml.dump(self.to_dict(), default_flow_style=None)
       
    def to_latex_doc(self):
        '''
        Parse the fields that need to be printed into latex.
        '''
        ### Latex
        latex_begin = '''\\documentclass[8pt]{report}
\\usepackage[a4paper,margin=0.1in,landscape]{geometry}
\\usepackage{amsmath}
\\usepackage{graphicx}

\\begin{document}
\\resizebox{0.98\\linewidth}{!}{%
'''
        latex_end = '''
}
\\end{document}
'''

        _eq_wrapper = lambda text: f'\\begin{{math}}\n{text}\n\\end{{math}}'
        latex = sympy.printing.latex

        ret_latex = latex_begin
        ret_latex += _eq_wrapper(latex(self.matrix))
        ret_latex += latex_end    

        return ret_latex
    

    def save_latex_doc(self, latex_filename: str):
        '''
        Save the Hamiltonian, its subspaces and eigenvectors in latex document
        '''

        ret_latex = self.parse_latex()
        with open(latex_filename, 'w') as ff:
            ff.write(ret_latex)

        return

    def _eq_wrapper(text):
        ret = f'\\begin{{math}}\n{text}\n\\end{{math}}'
        # if resizebox:
        #     ret = '\\resizebox{0.98\\linewidth}{!}{%\n'+ret+'\n}'

        return ret

    ##############################################################################
    # Numerical evaluation
    def _evaluate_str(self, formula: str, symbol_values: Dict[str, float]) -> float:
        '''
        Helper function to evaluate the numerical value 
        of a given formula
        '''

        translate_table = symbol_values
        translate_table.update(self._translate_table)

        return eval(formula, translate_table)
    
    def Neval_eigenvalues(self, symbol_values: Dict[str, float], shift_to_zero: bool=True, precision: int=10) -> List[Tuple[float, int]]:
        '''
        Numerically evaluate eigenvalues based on the `symbol values`
        substitution table. Will cast the eigenvalue to real numbers.

        Reises
        ------

        ValueError
            When the imaginary part of the eigenvalue is higher than 10^(-precision)
        '''
        res = []
        for formula, degeneracy in self.eigenvalues:
            value = np.around(self._evaluate_str(str(formula), symbol_values), precision)

            if np.abs(np.imag(value)) > 0:
                raise ValueError(f'Imaginary frequency {value}')
            
            res.append( (np.real(value), degeneracy) )

        if shift_to_zero:
            minE = min([E for E, _ in res])
            res = [(E-minE, deg) for E, deg in res]

        return res
    
    def Neval_eigenvectors(self, symbol_values: Dict[str, float], normalize: bool=False, precision: bool=10) -> List[List[float]]:
        '''
        Numerically evaluate eigenvectors based on the `symbol values`
        substitution table. Will cast the eigenvalue to real numbers.

        Reises
        ------

        ValueError
            When the imaginary part of the eigenvalue is higher than 10^(-precision)
        '''
        res = []
        for eigenspace in self.eigenvectors:
            eigenspace_num = []
            for eigenvector in eigenspace:
                ev_num = np.array([self._evaluate_str(str(x), symbol_values) for x in eigenvector], dtype=np.complex128)
                
                # Trim values
                ev_num = np.around(ev_num, precision)
                # Get rid of imag if smaller than precision
                if np.all( np.abs(np.imag(ev_num)) < np.power(1.0, -precision)):
                    ev_num = np.real(ev_num)

                if normalize:
                    ev_num /= np.linalg.norm(ev_num)

                eigenspace_num.append(ev_num)

            res.append(eigenspace_num)

        return res
    
    ##############################################################################
    # Diagonalization
    def determine_eigenvalues(self):
        '''
        Determine the eigenvalues of the Hamiltonian by constructing the
        characteristic polynomial and finding its roots.
        '''
        if self.eigenvalues:
            warnings.warn(f'Eigenvalues were already calculated. Redoing calculation.', UserWarning)
            self.eigenvalues = []    

        for eval_formula, eval_deg in sympy.roots(self.matrix.charpoly()).items():
            self.eigenvalues.append((eval_formula, eval_deg))

        return self.eigenvalues
    
    def determine_eigenvectors(self):
        '''
        Determine the eigenvectors of the Hamiltonian by determining the
        null space of the `H-eigenvalue*I` for each eigenvalue.
        '''
        if self.eigenvectors:
            warnings.warn(f'Eigenvectors were already calculated. Repeating calculation.', UserWarning)
            self.eigenvectors = []

        if not self.eigenvalues:
            warnings.warn(f'Eigenvalues not calculated in advance. Calculating eigenvalues first.')
            self.determine_eigenvalues()

        for eigenvalue_origin, (eigenvalue, degeneracy) in enumerate(self.eigenvalues):
            eigenvectors = (self.matrix - eigenvalue*sympy.eye(self.matrix.rank())).nullspace()

            for evec in eigenvectors:
                self.eigenvectors.append((evec, eigenvalue_origin))

        # for eval_formula, eval_deg in sympy.roots(self.matrix.charpoly()).items():
        #     self.eigenvalues.append((eval_formula, eval_deg))

        return self.eigenvectors

    ##############################################################################
    # Ordering the base
    def _sort_base(self):
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
        self._reorder_base(sorted_order)


    def _reorder_base(self, new_order: List[int]) -> Matrix:
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
    eigenvalues: list[...]
        List of eigenvalues with formula, degeneracy and subspace origin index.
    eigenvectors: list[...]
        List of eigenvectors ...
    '''
    def __init__(self, symmetry: str, Jval: float):
        self.Jval = Jval
        # Prepare a list of kets that form the basis
        if Jval%1==0:
            freeionkets = [f'{int(x)}' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]
        elif Jval%1==0.5:
            freeionkets = [f'{int(2*x):d}/2' for x in np.linspace(Jval,-Jval,int(2*Jval+1))]


        # Symmetry allowed Bnm according to McPhase manual
        # https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node133.html
        matrix_representation = self.H_from_point_group(symmetry, Jval)
        self.symmetry =  str(point_group_synonyms(symmetry)[0]) # + '\t(Hermann-Mauguin, Schoenflies, lattice names, id)'

        # Constructor of the parent class
        super().__init__(base=freeionkets, matrix=matrix_representation)

        # Other class specific fields
        self.subs = []

    def H_from_point_group(self, point_group_name: str, Jval: float) -> sympy.Matrix:
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
    def determine_eigenvalues(self, from_blocks: bool=True):
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
    
    def determine_eigenvectors(self, from_blocks: bool=True):
        '''
        Updated fields
        --------------
        subs
        '''
        if not from_blocks:
            return super().determine_eigenvectors()
        
        if not self.subs:
            self.make_block_form()

        if not self.eigenvalues:
            self.determine_eigenvalues()

        # This approach looses information about the subspace origin
        # and does not make eigenvectors from complete Hamiltonian 
        # Hilbert space
        self.eigenvectors = []
        for Hsub in self.subs:
            Hsub.determine_eigenvectors()
            self.eigenvectors.append(Hsub.eigenvectors)

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

    def to_dict(self):
        '''
        Represent as dictionary
        '''

        H_dict = {}

        H_dict['base'] = self.base
        H_dict['J'] = self.Jval
        H_dict['symmetry'] = self.symmetry
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

        # Eigenvectors -> ???
        # H_dict['eigenvectors'] = {}
        # for n, (evec, eval_origin) in enumerate(self.eigenvectors):
        #     eigenvector = [str(u) for u in evec]
        #     H_dict['eigenvectors'][n] = dict(formula=eigenvector, eval_origin=eval_origin)

        H_dict['subspaces'] = {}
        for n, Hsub in enumerate(self.subs):
            H_dict['subspaces'][n] = Hsub.to_dict()

        return H_dict

    # def to_yaml(self):
    #     '''
    #     Parse the fields to YAML
    #     '''

    #     ret = '---\n'

    #     # Base
    #     ret += 'base: [' + ', '.join(self.base) + ']\n'

    #     # Hamiltonian
    #     ret += 'hamiltonian:\n'
    #     hamiltonian_str = str(self.matrix)
    #     it_i = hamiltonian_str.find('[[') + 2
    #     it_f = hamiltonian_str.find(']]')
    #     for row in hamiltonian_str[it_i: it_f].split('], ['):
    #         ret += '  - [' + str(row) + ']\n'

    #     # Blocks
    #     for Hsub in self.subs:
    #         ret += '\n'
    #         ret += Hsub.to_yaml()

    #     # # Eigenvalues
    #     # ret += 'eigenvalues:'
    #     # if self.eigenvalues:
    #     #     ret += '\n'
    #     #     for eval in self.eigenvalues:
    #     #         ret += '  - ' + str(eval) + '\n'
    #     # else:
    #     #     ret += ' []'

    #     return ret

    def to_latex_doc(self):
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
    symmetry = 'C2v'

    H = CEF_Hamiltonian(symmetry=symmetry, Jval=Jval)
    # print(H)
    # print(H.to_yaml_auto())
    # print(H.to_yaml())

    H.make_block_form()
    # print(H.to_yaml_auto())
    # print(H.to_yaml())

    evals = H.determine_eigenvalues()
    print(H.eigenvalues)
    print(H.to_yaml_auto())
    # print(H.to_yaml())

    # evals = H.determine_eigenvectors()
    # print(H.to_yaml_auto())
    # print(H.to_yaml())

    # print('Loaded:')
    # data = yaml.load(H.to_yaml(), Loader=yaml.Loader)
    # for t in data:
    #     print(t, data[t])

    # print('Dicted:')
    # print(H.to_dict())