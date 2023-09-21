    ##############################################################################33
    # Printing and representing
    def parse_latex(self):
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
    

    def save_latex(self, latex_filename: str):
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
