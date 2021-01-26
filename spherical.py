import numpy as np
from scipy.special import comb, factorial, factorial2

def shift(m):

    if m < 0:
        return 0.5
    else:
        return 0.0

def shCoefficient(l, m, i, j, k):
    #coefficients of spherical harmonics 
    factor = pow(0.25,i) * comb(l, i) * comb(l - i, abs(m) + i) * comb(i, j) * comb(abs(m), 2 * k)

    if m < 0:
        return np.real(pow(complex(-1), i + k - shift(m))) * factor
    else:
        return pow(-1, i + k - shift(m)) * factor
    

def shNormal(l, m):
 
    return (1 / (pow(2, abs(m)) * factorial(l))) * np.sqrt((2 * factorial(l+ abs(m)) * factorial(l - abs(m)))
        / pow(2 ,int(m == 0))    )


def sh(l, m):

    spherical = {} 
    normal = shNormal(l, m)

    for i, j, k in [ (x, y, z + shift(m)) for x in range(((l - abs(m)) // 2) + 1)
        for y in range(x + 1) for z in range((abs(m) // 2) + 1)  ]:

        coefficient = shCoefficient(l, m, i, j, k) * normal
        # a_x, a_y, a_z are the components of a Cartesian Gaussian primitive
        a_x = int(2 * i + abs(m) - (2 * (j + k)))
        a_y = int(2 * (j + k))
        a_z = int(l - 2 * i - abs(m))
        if (coefficient != 0 ) and (coefficient in spherical):
            spherical[(a_x, a_y, a_z)] += coefficient
        elif coefficient != 0 :
            spherical[(a_x, a_y, a_z)] = coefficient



    return dict(spherical)

def shTransformation(l, cartesianOrder, sphericalOrder, direction):


    order = {
             components: index for index, components in enumerate(list(map(tuple, cartesianOrder)))
            }
    transform = np.zeros(((l + 1) * (l + 2) // 2, 2 * l + 1))

    for i, shFunction in enumerate(sphericalOrder):
        if shFunction[0] == "-":
            sign = -1
            shFunction = shFunction[1:]
        else:
            sign = 1
        if shFunction[0] == "s":
            mag = -1 * int(shFunction[1:])
        else:
            mag = int(shFunction[1:])
        harmonic = sh(l, mag)
        for components, coeff in harmonic.items():
            transform[order[components], i] = sign * coeff

    transform *= (np.sqrt(np.prod(factorial2(2 * cartesianOrder - 1), axis=1))[:, np.newaxis]) / np.sqrt(factorial2(2 * l - 1))

    if direction == "c->s":
        return transform.T
    return transform


print(shTransformation(2,np.asarray([[2,0,0],[1,1,0],[1,0,1],[0,2,0],[0,1,1],[0,0,2]]),['c0','c1','s1','c2','s2'],'c->s'))


#==========================================Hyperpolarizability
import numpy as np
import rhf
from post import polarizabilities
from basis import electronCount

molAtom, molBasis, molData = rhf.mol([])
energy = rhf.scf(molAtom, molBasis, molData,[])

data = np.load('../mints/h2o-sto-3g-mints.npz')
c = data['c']
f = data['f']
ERI = data['i']
charge = molData[0]
nOccupied = electronCount(molAtom, charge) //2
j = data['j']
#---------------------------------------------------------

def hyperPolarizabilities(molAtom, molBasis, c, f, j, ERI, nOccupied):

	#get responses and dipole tensors
	responses , dipoles_mo = polarizabilities(molAtom, molBasis, c, f, ERI, nOccupied)

	#number of components
	nComponents = responses.shape[0]

	nOrbitals = c.shape[0]
	nVirtual = nOrbitals - nOccupied

	#reshape responses (3, vo) -> (3,o+v,o+v) compatible with dipoles
	u = np.zeros_like(dipoles_mo)
	uv = np.zeros((nComponents, nOccupied, nVirtual))
	for comp in range(nComponents):
		k = 0
		for p in range(nOccupied):
			for q in range(nVirtual):
				u[comp, p, nOccupied+q] = 0.5 * responses[comp, k]
				uv[comp, p, q] = responses[comp, k]
				u[comp, q + nOccupied, p] = -0.5 * responses[comp, k]
				k += 1

	#get occupied and virtual c matrices
	co = np.zeros((nOrbitals, nOccupied))
	co = c[:, :nOccupied]
	cv = np.zeros((nOrbitals, nVirtual))
	cv = c[:,nOccupied:]
	
	r = np.zeros((nOrbitals, nOccupied))

	for comp in range(nComponents):

		r = np.dot(uv[comp], cv.T).T
			
		d = 2 * np.dot(co, r.T)
		print(np.dot(d, j))


hyperPolarizabilities(molAtom, molBasis, c, f, j, ERI, nOccupied)
