# Tools
import numpy as np
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
from .data_containers import sympy_expression_mapping


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
        Lsit of eigenvalues, where each entry contains eigenvalue's formula and degeneracy.
    eigenvectors: List[Tuple[expr, int]]
        List of eigenvectors, where each entry contains the eigenvector's formula
        and index of the eigenvalue from which it originates.

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
    
    def to_dict(self) -> Dict[str, Any]:
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
    
    def to_yaml(self) -> str:
        '''
        Convert the Hamiltonian to YAML format.

        Works by calling the `self.to_dict()` function.
        '''
        return yaml.dump(self.to_dict(), default_flow_style=None)
       
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

    def _eq_wrapper(text: str) -> str:
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
    def determine_eigenvalues(self) -> List[Tuple[expr, int]]:
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
    
    def determine_eigenvectors(self) -> List[Tuple[expr, int]]:
        '''
        Determine the eigenvectors of the Hamiltonian by determining the
        null space of the `H-eigenvalue*I` for each eigenvalue.

        Returns
        -------
            List of tuples, where each tuple contains the formular of the eigenvector
            and the index of the eigenvalue to which it corresponds.
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