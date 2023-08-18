'''
Most important things so far:
 1. The eigenvector convention is that the first entry corresponds to highest spin, i.e.
    for the J=3/2 system |3/2> = [1,0,0,0].
 2. The CEF operators correspond to the Stevens notation.
 3. Stevens operators correspond to observables, thus they should be hermitian.
    By running the library a check is performed to see if all Onm matrices are hermitian.
'''

from sympy import Integer, Rational, I, sqrt, eye, diag, conjugate, transpose, zeros


def J_z(J):
	return diag(*[Rational(2*(J-m),2) for m in range(Integer(2*J)+1)])

def J_y(J):
	return -I * Rational(1,2) * (J_plus(J) - J_minus(J))

def J_x(J):
	return Rational(1,2) * (J_plus(J) + J_minus(J))


def J_plus(J):
	J = Rational(2*J,2)
	J_plus = zeros(2*J+1, 2*J+1)
	for n in range(2*J):
		m = J-n-1
		J_plus[n, n+1] = sqrt(J*(J+1) - m*(m+1))

	return J_plus

def J_minus(J):
	return transpose(conjugate(J_plus(J)))



def O_20(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	Jz = J_z(J)
	E = eye(J2p1)

	return Integer(3)*Jz*Jz - JJ * E

def O_22(J):
	Jplus = J_plus(J)
	Jminus = J_minus(J)

	return Rational(1,2) * (Jplus**2 + Jminus**2)

def O_2m2(J):
	Jplus = J_plus(J)
	Jminus = J_minus(J)

	return -I * Rational(1,2) * (Jplus**2 - Jminus**2)

def O_40(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	Jz = J_z(J)
	E = eye(J2p1)

	return 35*Jz**4 + (25 - 30 * JJ)*Jz**2 + E*(3 * JJ**2 - 6 * JJ)


def O_42(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)
	E = eye(J2p1)

    # helper matrices:	
	M_1 = 7 * Jz**2 - E*(JJ + 5)
	M_2 = Jplus**2 + Jminus**2

	return Rational(1,4) * (M_1@M_2 + M_2@M_1)

def O_4m2(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)
	E = eye(J2p1)

    # helper matrices:	
	M_1 = 7*Jz**2 - E*(JJ + 5)
	M_2 = Jplus**2 - Jminus**2

	return -I * Rational(1,4) * (M_1@M_2 + M_2@M_1)

def O_43(J):
	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

    # helper matrices:	
	M_1 = Jz
	M_2 = Jplus**3 + Jminus**3
	
	return Rational(1,4) * (M_1@M_2 + M_2@M_1)

def O_4m3(J):
	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

    # helper matrices:	
	M_1 = Jz
	M_2 = Jplus**3 - Jminus**3
	
	return -I * Rational(1,4) * (M_1@M_2 + M_2@M_1)

def O_44(J):
	Jplus = J_plus(J)
	Jminus = J_minus(J)
	
	return Rational(1,2) * (Jplus**4 + Jminus**4)

def O_4m4(J):
	Jplus = J_plus(J)
	Jminus = J_minus(J)
	
	return -I*Rational(1,2) * (Jplus**4 - Jminus**4)

def O_60(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	Jz = J_z(J)
	E = eye(J2p1)

	return 231*Jz**6 + Jz**4 * (735 - 315 * JJ) + Jz**2*(105*JJ**2 - 525*JJ + 294) + E*(-5*JJ**3 + 40 * JJ**2 - 60*JJ)

def O_62(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	E = eye(J2p1)
	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

	# helper matrices:	
	M_1 = 33*Jz**4 - Jz**2 * (18 * JJ + 123) + E*(JJ**2 + 10*JJ + 102)
	M_2 = Jplus**2 + Jminus**2

	return Rational(1,4) * (M_1 @ M_2 + M_2 @ M_1)

def O_6m2(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	E = eye(J2p1)
	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

	# helper matrices:	
	M_1 = 33*Jz**4 - Jz**2 * (18 * JJ + 123) + E*(JJ**2 + 10*JJ + 102)
	M_2 = Jplus**2 - Jminus**2

	return -I*Rational(1,4) * (M_1 @ M_2 + M_2 @ M_1)

def O_63(J):
	JJ = Rational(2*J*(J+1), 2)

	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

	# helper matrices:	
	M_1 = 11*Jz**3 - Jz*(59 + 3*JJ)
	M_2 = Jplus**3 + Jminus**3

	return Rational(1,4) * (M_1 @ M_2 + M_2 @ M_1)

def O_6m3(J):
	JJ = Rational(2*J*(J+1), 2)

	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

	# helper matrices:	
	M_1 = 11*Jz**3 - Jz*(59 + 3*JJ)
	M_2 = Jplus**3 - Jminus**3

	return -I*Rational(1,4) * (M_1 @ M_2 + M_2 @ M_1)

def O_64(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	E = eye(J2p1)
	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

    	# helper matrices:	
	M_1 = 11*Jz**2 - E*(JJ + 38)
	M_2 = Jplus**4 + Jminus**4

	return Rational(1,4) * (M_1 @ M_2 + M_2 @ M_1)

def O_6m4(J):
	JJ = Rational(2*J*(J+1), 2)
	J2p1 = Integer(2*J + 1)

	E = eye(J2p1)
	Jz = J_z(J)
	Jplus = J_plus(J)
	Jminus = J_minus(J)

    	# helper matrices:	
	M_1 = 11*Jz**2 - E*(JJ + 38)
	M_2 = Jplus**4 - Jminus**4

	return -I*Rational(1,4) * (M_1 @ M_2 + M_2 @ M_1)


def O_66(J):
	Jplus = J_plus(J)
	Jminus = J_minus(J)
	
	return Rational(1,2) * (Jplus**6 + Jminus**6)

def O_6m6(J):
	Jplus = J_plus(J)
	Jminus = J_minus(J)
	
	return -I*Rational(1,2) * (Jplus**6 - Jminus**6)


if __name__ == '__main__':
	import cef_matrices
	from sympy import pprint
	Jval = 4.5

	for func_name in dir(cef_matrices):
		if func_name[0]=='O' or func_name=='J_x' or func_name=='J_y':
			print(f"Testing matrix {func_name} for J={Jval}")
			matrix_function = getattr(cef_matrices, func_name)
			matrix = matrix_function(Jval)

			pprint(matrix)

			if not matrix==transpose(conjugate(matrix)):
				raise ValueError(f'Matrix {func_name} is not hermitian')