import sys
PATH = 'C:/Users/Stekiel/Documents/GitHub/crysfipy_symbolic'
sys.path.append(PATH)

from crysfipy_symbolic.cefhamiltonian import CEF_Hamiltonian
from crysfipy_symbolic.data_containers import point_group_synonyms

Jval = 2.0
lattice_type = 'cubic_1'

print(f'Constructing tables for J={Jval}, {lattice_type}')

H = CEF_Hamiltonian(symmetry=lattice_type, Jval=Jval)
H.make_block_form()
H.determine_eigenvectors()

print('Finished diagonalization...')

Jstr = str(int(Jval*10)/10).replace('.', 'p')
table_name = f'{lattice_type}__J_{Jstr}'

# print('Saving TEX file...')
# with open(f'{PATH}/tables/{table_name}.tex', 'w') as ff:
#     ff.write(H.to_latex_doc())

print('Saving YAML file...')
with open(f'{PATH}/tables/{table_name}.yaml', 'w') as ff:
    ff.write(H.to_yaml())