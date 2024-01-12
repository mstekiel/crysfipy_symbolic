import numpy as np
from typing import Union

def sympy_expression_mapping():
    return {'sqrt':np.emath.sqrt, 'I':1j}


def point_group_synonyms(name: Union[str, int]):
    '''
    Return the synonyms of a given point group name.

    Contains database of Hermann-Mauguin, Schoenflies and lattice names as np.array with
    dtype = [
            ('HM', 'U8'),
            ('Sch', 'U8'),
            ('lattice', 'U16'),
            ('id', 'i4')
            ]

    Parameters:
    -----------
    name:
        Name of the space group in Hermann-Mauguin, Schoenflies or lattice type notation.
        Alternatively, integer identifyer of the point symmetry.

        
    Returns:
    --------
    Tuple of (Hermann-Mauguin, Schoenflies, lattice names, id),
    in case the name is ambiguous a list of fitting tuples.
    '''


    symmetry_synonyms = [
        # '0' reserved for future triclinic
        ('2',   'C2',   'monoclinic',   1),
        ('m',   'Cs',   'monoclinic',   1),
        ('2/m', 'C2h',  'monoclinic',   1),
        ('mm2', 'C2v',  'orthorhombic', 2),
        ('222', 'D2',   'orthorhombic', 2),
        ('mmm', 'D2h',  'orthorhombic', 2),
        ('4',   'C4',   'tetragonal_1', 3),
        ('-4',  'S4',   'tetragonal_1', 3),
        ('4/m', 'C4h',  'tetragonal_1', 3),
        ('422', 'D4',   'tetragonal_2', 4),
        ('4mm', 'C4v',  'tetragonal_2', 4),
        ('-42m',        'D2d',  'tetragonal_2', 4),
        ('4/mmm',       'D4h',  'tetragonal_2', 4),
        ('3',   'C3',   'trigonal_1',   5),
        ('-3',  'S6',   'trigonal_1',   5),
        ('32',  'D3',   'trigonal_2',   6),
        ('3m',  'C3v',  'trigonal_2',   6),
        ('-3m', 'D3d',  'trigonal_2',   6),
        ('6',   'C6',   'hexagonal_1',  7),
        ('-6',  'C3h',  'hexagonal_1',  7),
        ('6/m', 'C6h',  'hexagonal_1',  7),
        ('622', 'D6',   'hexagonal_2',  8),
        ('6mm', 'C6v',  'hexagonal_2',  8),
        ('-6m2',        'D3h',  'hexagonal_2',  8),
        ('6/mmm',       'D6h',  'hexagonal_2',  8),
        ('23',  'T',    'cubic_1',      9),
        ('m-3', 'Th',   'cubic_1',      9),
        ('-43m',        'Td',   'cubic_2',      10),
        ('432', 'O',    'cubic_2',      10),
        ('m-3m',        'Oh',   'cubic_2',      10)
    ]

    symmetry_db = np.array(symmetry_synonyms, dtype=[
                                        ('HM', 'U8'),
                                        ('Sch', 'U8'),
                                        ('lattice', 'U16'),
                                        ('id', 'i4')
                                        ])
    
    symmetry = None

    if isinstance(name, str):
        if name in symmetry_db['HM']:
            symmetry = symmetry_db[ symmetry_db['HM']==name ]
        if name in symmetry_db['Sch']:
            symmetry = symmetry_db[ symmetry_db['Sch']==name ]
        if name in symmetry_db['lattice']:
            symmetry = symmetry_db[ symmetry_db['lattice']==name ]
    if isinstance(name, int):
        if name in symmetry_db['id']:
            symmetry = symmetry_db[ symmetry_db['id']==name ]

    if symmetry is None:
        raise NameError(f'{name} is an unknown point group identifier.')

    return symmetry

# class CEF_tables():
# 	pass
# 	'''
# 	Hard-coded data on various CEF related objects, printing and conversion functionalities.
#
# 	'''
# 	def __init__(self):
# 		return

# 	def fill_symmetry(self):
# 		return


if __name__ == '__main__':
    db = point_group_synonyms(6)

    print(db)