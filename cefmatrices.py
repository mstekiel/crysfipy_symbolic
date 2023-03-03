'''
Additions by Michal Stekiel

The matrices and conventions need to be cross-checked against some reliable book.
Most important things so far:
 1. The eigenvector convention is that the first entry corresponds to highest spin, i.e.
    for the J=3/2 system |3/2> = [1,0,0,0].
    The convention can be flipped so that first entry corresponds to the lowest spin,
    by the 'convention=-1' keyword in all matrices.
 2. The CEF operators correspond to the Stevens notation.
'''

import sympy
from numpy import diag, linspace, conj, transpose, sqrt, eye, dot, array
from numpy.linalg import matrix_power as mp


def J_z(J,convention=1):
	return diag([convention*sympy.Rational(2*(J-m),2) for m in range(2*J+1)])
    # return diag(linspace(convention * J,convention * (-J),int(2*J+1)))

def J_y(J,convention = 1):
	return .5/1.j*(J_plus(J,convention) - J_minus(J,convention))

def J_x(J,convention = 1):
	return .5 * (J_plus(J,convention) + J_minus(J,convention))



def J_plus(J,convention = 1):
	p1 = [sympy.sqrt(J*(J+1) + sympy.Rational(2*(J-m),2)*(1-sympy.Rational(2*(J-m),2))) for m in range(2*J)]

	# p1 = linspace(-J,J-1,int(2*J))
	# p1 = sqrt(J*(J+1) - p1*(p1+1))

	return diag(p1,convention)

def J_minus(J,convention = 1):
	return J_plus(J,convention).conj().transpose()



def O_20(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jz = J_z(J,convention)
	E = eye(J2p1)

	return 3*Jz*Jz - JJ * E

def O_22(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)

	return 0.5 * (mp(Jplus, 2) + mp(Jminus, 2))

def O_40(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jz = J_z(J,convention)
	E = eye(J2p1)

	return 35 * mp(Jz, 4) + (25 - 30 * JJ)*mp(Jz, 2) + E * JJ * (3 * JJ - 6)


def O_42(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jz = J_z(J,convention)
	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)
	E = eye(J2p1)

    # helper matrices:	
	M_1 = 7 * mp(Jz, 2) - E*(JJ + 5)
	M_2 = mp(Jplus, 2) + mp(Jminus, 2)

	return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))

def O_43(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jz = J_z(J,convention)
	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)

    # helper matrices:	
	M_1 = Jz
	M_2 = mp(Jplus, 3) + mp(Jminus, 3)
	
	return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))


def O_44(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)
	
	return 0.5 * (mp(Jplus, 4) + mp(Jminus, 4))


def O_60(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jz = J_z(J,convention)
	E = eye(J2p1)

	return 231 * mp(Jz, 6) + mp(Jz, 4) * (735 - 315 * JJ) + mp(Jz, 2)*(105*JJ**2 - 525*JJ + 294) + dot(E, -5*JJ**3 + 40 * JJ**2 - 60*JJ)

def O_62(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)
	
	E = eye(J2p1)
	Jz = J_z(J,convention)
	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)
	
    # helper matrices:	
	M_1 = 33 * mp(Jz, 4) - mp(Jz, 2) * (18 * JJ + 123) + dot(E, JJ**2 + 10*JJ + 102)
	M_2 = mp(Jplus, 2) + mp(Jminus, 2)

	return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))

def O_63(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jz = J_z(J,convention)
	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)

    # helper matrices:	
	M_1 = 11 * mp(Jz, 3) - Jz * (59 + 3*JJ)
	M_2 = mp(Jplus, 3) + mp(Jminus, 3)
	
	return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))

def O_64(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	E = eye(J2p1)
	Jz = J_z(J,convention)
	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)

    # helper matrices:	
	M_1 = 11 * mp(Jz, 2) - E * (JJ + 38)
	M_2 = mp(Jplus, 4) + mp(Jminus, 4)

	return 0.25 * (dot(M_1, M_2) + dot(M_2, M_1))


def O_66(J,convention = 1):
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jplus = J_plus(J,convention)
	Jminus = J_minus(J,convention)
	
	return 0.5 * (mp(Jplus, 6) + mp(Jminus, 6))


if __name__ == '__main__':
	J = 2
	JJ = J*(J+1)
	J2p1 = int(2*J + 1)

	Jz = J_z(J,1)
	E = eye(J2p1)

	O_60_test =  231 * mp(Jz, 6) + mp(Jz, 4) * (735 - 315 * JJ) + mp(Jz, 2)*(105*JJ**2 - 525*JJ + 294) + dot(E, -5*JJ**3 + 40 * JJ**2 - 60*JJ)

	print(231 * mp(Jz, 6))
	print(mp(Jz, 4) * (735 - 315 * JJ))
	print(mp(Jz, 2)*(105*JJ**2 - 525*JJ + 294))
	print(dot(E, -5*JJ**3 + 40 * JJ**2 - 60*JJ))
	print(O_60_test)