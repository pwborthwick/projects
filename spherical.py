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


