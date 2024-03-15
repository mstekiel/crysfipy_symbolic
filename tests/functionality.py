import sys
PATH = 'C:/Users/Stekiel/Documents/GitHub/crysfipy_symbolic'
sys.path.append(PATH)

from crysfipy_symbolic.cefhamiltonian import CEF_Hamiltonian
from crysfipy_symbolic.data_containers import point_group_synonyms


import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

symmetry = 'mmm'
Jval=2

H = CEF_Hamiltonian(symmetry=symmetry, Jval=Jval)
H.make_block_form()
H.determine_eigenvectors()

print(H.to_yaml())