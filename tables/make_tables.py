import sys
PATH = 'C:/Users/Stekiel/Documents/GitHub/crysfipy_symbolic'
sys.path.append(PATH)

from crysfipy_symbolic import CEF_Hamiltonian
from crysfipy_symbolic import point_group_synonyms

import numpy as np
from timeit import default_timer as timer

H = CEF_Hamiltonian(symmetry=1, Jval=2.5)


symmetry_range = range(2, 10)
Jval_range = np.arange(1, 10 + 0.1, 0.5)
t_cutoff = 8

results = []

### With eigenvalues
for symmetry in symmetry_range:
    t_sofar = []
    for Jval in Jval_range:
        HM, Sch, lattice_type, _ = point_group_synonyms(symmetry)[0]
        Jstr = str(int(Jval*10)/10).replace('.', 'p')

        table_name = f'{lattice_type}__J_{Jstr}'
        print(f'J={Jval}, {lattice_type}, +eigenvectors')

        start_time = timer()
        H = CEF_Hamiltonian(symmetry=HM, Jval=Jval)
        H.make_block_form()
        try:
            H.determine_eigenvectors()

            print('Saving YAML file...')
            with open(f'{PATH}/tables/{table_name}.yaml', 'w') as ff:
                ff.write(H.to_yaml())
        except ValueError:
            print(f'WARNING: J={Jval} {lattice_type} Hamiltonian zero')
            continue
        finish_time = timer()

        calculation_time = finish_time - start_time

        print(f'\t t = {calculation_time:.3} s')
        results.append([Jval, lattice_type, calculation_time, True])
        t_sofar.append(calculation_time)

        t_next = 1.0
        if len(t_sofar) > 1:
            acceleration = np.mean([tn/tp for tp,tn in zip(t_sofar[:-1], t_sofar[1:])])
            t_next = calculation_time*acceleration

        print(f'\t t_next = {t_next:.3} s')
        if t_next > t_cutoff:
            break
