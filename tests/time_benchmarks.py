import sys
PATH = 'C:/Users/Stekiel/Documents/GitHub/crysfipy_symbolic'
sys.path.append(PATH)

from crysfipy_symbolic.cefhamiltonian import CEF_Hamiltonian
from crysfipy_symbolic.data_containers import point_group_synonyms


import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer



symmetry_range = range(2, 10)
Jval_range = np.arange(1, 10 + 0.1, 0.5)
t_cutoff = 1

results = []

#########################################################################
### Without eigenvalues
for symmetry in symmetry_range:
    t_sofar = []
    for Jval in Jval_range:
        HM, Sch, lattice_type, _ = point_group_synonyms(symmetry)[0]
        print(f'J={Jval}, {lattice_type}')

        start_time = timer()
        H = CEF_Hamiltonian(symmetry=HM, Jval=Jval)
        H.make_block_form()
        try:
            H.determine_eigenvalues()
        except ValueError:
            print(f'WARNING: J={Jval} {lattice_type} Hamiltonian zero')
            continue
        finish_time = timer()

        calculation_time = finish_time - start_time

        print(f'\t t = {calculation_time:.3} s')
        results.append([Jval, lattice_type, calculation_time, False])
        t_sofar.append(calculation_time)

        t_next = calculation_time
        if len(t_sofar) > 2:
            acceleration = np.mean([tn/tp for tp,tn in zip(t_sofar[:-1], t_sofar[1:])])
            t_next = calculation_time*acceleration

        print(f'\t t_next = {t_next:.3} s')
        if t_next > t_cutoff:
            break

#########################################################################
### With eigenvalues
for symmetry in symmetry_range:
    t_sofar = []
    for Jval in Jval_range:
        HM, Sch, lattice_type, _ = point_group_synonyms(symmetry)[0]
        print(f'J={Jval}, {lattice_type}, +eigenvectors')

        start_time = timer()
        H = CEF_Hamiltonian(symmetry=HM, Jval=Jval)
        H.make_block_form()
        try:
            H.determine_eigenvectors()
        except ValueError:
            print(f'WARNING: J={Jval} {lattice_type} Hamiltonian zero')
            continue
        finish_time = timer()

        calculation_time = finish_time - start_time

        print(f'\t t = {calculation_time:.3} s')
        results.append([Jval, lattice_type, calculation_time, True])
        t_sofar.append(calculation_time)

        t_next = 1.0
        if len(t_sofar) > 2:
            acceleration = np.mean([tn/tp for tp,tn in zip(t_sofar[:-1], t_sofar[1:])])
            t_next = calculation_time*acceleration

        print(f'\t t_next = {t_next:.3} s')
        if t_next > t_cutoff:
            break

for result in results:
    print(result)

#########################################################################
### Plot
fig, ax = plt.subplots()

ax.set_xlabel('J')
ax.set_ylabel('diagonalization time (s)')
ax.set_yscale('log')

for lattice_type in np.unique([r[1] for r in results]):
    # Without eigenvectors
    J, time = np.transpose([[r[0], r[2]] for r in results if r[1]==lattice_type and not r[3]])
    line = ax.plot(J, time, 'o-', label=lattice_type)

    # With eigenvectors
    J, time = np.transpose([[r[0], r[2]] for r in results if r[1]==lattice_type and r[3]])
    ax.plot(J, time, 's:', color=line[0].get_color(), label=lattice_type+'+ev')

ax.legend()
fig.tight_layout()
fig.savefig(f'{PATH}/tests/diagonalization_time-short.png', dpi=100)
