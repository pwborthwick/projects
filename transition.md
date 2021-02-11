## Transition Dipoles

I want to get the transition dipoles, there is a psi4 jupyter notebook [here](https://hub.gke2.mybinder.org/user/robertodr-tddft_psicon2020-p9gjw54v/lab) which I shall base this on. Firstly 

### Dipoles
The example is H<sub>2</sub>O in the STO-3G basis. \
From mints get **c** eigenvectors, then generate **c**<sub>o</sub> as c[:, :n<sub>o</sub>] and **c**<sub>v</sub> as c[:, n<sub>o</sub>:], where n<sub>o</sub> is the number of (doubly) occupied orbitals. \
Writing n<sub>e</sub> is number of electrons, n<sub>b</sub> is the number of basis functions, n<sub>v</sub> is the number of virtual orbitals and n<sub>s</sub> is the number of spin orbitals (= 2n<sub>b</sub>). \
In harpy this is implemented as
```python
mints = np.load('../mints/h2o-sto-3g-mints.npz')

Ne, Nb, No, Nv, Ns = [10, 7, 5, 2, 20]

c = mints['c']
co = c[:, :No]
cv = c[:, No:]
```
The dipoles are a [n<sub>b</sub>,n<sub>b</sub>] matrix (**&#956;**), to get to the right shape form \
**c**<sub>o</sub><sup>T</sup>**.&#956;.c**<sub>v</sub>, this will have shape [n<sub>o</sub>, n<sub>b</sub>][n<sub>b</sub>, n<sub>b</sub>][n<sub>b</sub>, n<sub>v</sub>] -> [n<sub>o</sub>, n<sub>v</sub>] \
In harpy this is implemented as
```python
from post import dipoleComponent

atoms, bases, data = rhf.mol([])
mu = np.zeros((3, Nb, Nb))

axes = ['x', 'y', 'z']
for i, axis in enumerate(axes):
    mu[i] = dipoleComponent(atoms, bases, axis, 'nuclear charge')

muMO = np.zeros((3, No, Nv))
for i in range(3):
    muMO[i] = np.dot(co.T, np.dot(mu[i], cv))
```
From psi4 we get for mu[0] 
```
 [-2.34699852e-15 -5.47538223e-02]
 [ 1.11162392e-14 -7.75544038e-02]
 [-9.15705544e-01 -4.08772620e-14]
 [-5.47495359e-14  7.80910397e-01]
 [-5.28535409e-15 -1.05919482e-14]
```
and from harpy
```
 [-1.56125113e-16 -5.47538454e-02]
 [ 1.83186799e-15 -7.75543066e-02]
 [-9.15706252e-01  2.16493490e-15]
 [ 3.66373598e-15  7.80909455e-01]
 [ 1.48260384e-15  7.79000144e-30]
```
so we're OK here.

###alpha and beta eigenvectors
To get these we need to do a TDHF calculation. To avoid the left and right eigenvector problem for the moment we'll use the Hermitian form of the equations \
The basic form of the equations \
**A** = &#948;<sub>𝑖𝑗</sub>**𝑓**<sub>𝑎𝑏</sub> − &#948;<sub>𝑎𝑏</sub>**𝑓**<sub>𝑖𝑗</sub> + <bi||𝑗a> \
**B** = <bi||ja> \

(**A** - **B**)(**A** + **B**) |**X**<sub>n</sub> + **Y**<sub>n</sub>> =  	&#969;<sup>2</sup> |**X**<sub>n</sub> + **Y**<sub>n</sub>> (i)\
left multiply by (**A** - **B**)<sup>-&#189;</sup> and right multiply by (**A** - **B**)<sup>&#189;</sup> gives \
(**A** - **B**)<sup>&#189;</sup>**.**(**A** + **B**)**.**(**A** - **B**)<sup>&#189;</sup> **T**<sub>n</sub> = &#969;<sup>2</sup> **T**<sub>n</sub> (ii) \
writing (ii) as (left multiply by (**A** - **B**)<sup>&#189;</sup>) gives \
(**A** - **B**)(**A** + **B**) ((**A** - **B**)<sup>&#189;</sup> 
**T**<sub>n</sub>) =  &#969;<sup>2</sup> ((**A** - **B**)<sup>&#189;</sup> **T**<sub>n</sub>) \
comparing this Hermitian form with (i) gives  |**X**<sub>n</sub> + **Y**<sub>n</sub>>  =  ((**A** - **B**)<sup>&#189;</sup> **T**<sub>n</sub>) (iii) \

Running a TDHF calculation using the Hermitian form gives eigenvalues \
```
[  7.75970415   7.75970415   7.75970415   8.15643883   8.15643883
   8.15643883   9.59546498   9.59546498   9.59546498   9.65401306
   9.93573411   9.93573411   9.93573411  11.30137037  13.60845200 ]
```
from psi4
```
[9.654017765570021  11.301372582481267  13.608453093245776]
```
This for first 3 singlets - OK so far.

The eigenvectors from harpy psi4 for alpha eigenvectors right for the first and second singles singlets - alpha and beta are the same
```
[[ 1.24887231e-16  5.76788541e-19]
 [-1.08892241e-15 -6.05807431e-17]
 [ 2.13542747e-15 -5.54872337e-16]
 [ 3.94632694e-15  3.10090564e-16]
 [ 9.52476119e-01  1.79097524e-15]]
 
[[ 2.04784157e-16  1.91983878e-17]
 [ 2.18676394e-14 -4.43304138e-16]
 [-3.57368101e-14  1.58127427e-15]
 [ 1.34315726e-13 -1.89111781e-15]
 [-2.64227345e-14  9.70320358e-01]]
```
The harpy implementstion is 
```
spinOrbitals = len(bases) * 2
#get fock matrix in MO spin basis
fockMOspin = buildFockMOspin(spinOrbitals, eigenVectors, fock)
  
#get two-electron repulsion integrals in MO basis
eriMO = buildEriMO(eigenVectors, ERI)

#transform eri from MO to spin basis
eriMOspin = buildEriDoubleBar(spinOrbitals, eriMO)

#orbital occupation
nElectrons = Ne
nOccupied  = No
nVirtual   = Nv

A = np.zeros((nOccupied*nVirtual*4, nOccupied*nVirtual*4))
B = np.zeros((nOccupied*nVirtual*4, nOccupied*nVirtual*4)) 

ia = 0
for i in range(0, nElectrons):
  for a in range(nElectrons, spinOrbitals):

    jb = 0
    for j in range(0, nElectrons):
      for b in range(nElectrons, spinOrbitals):

        if (i == j):
          A[ia,jb] += fockMOspin[a,b]
        if (a == b):
          A[ia,jb] -= fockMOspin[i,j]

        A[ia,jb] += eriMOspin[a,j,i,b]
        B[ia,jb] += eriMOspin[a,b,i,j]

        jb += 1
    ia += 1

HA = fractPow(A - B, 0.5)
HB = A + B
rpaH = np.dot(HA,np.dot(HB, HA))

from scipy.linalg import lu, eig
w , vr = np.linalg.eigh(rpaH)
#solve T for |X+Y>
vr = np.dot(HA, vr)            

print('\neigenvalues...')
print(np.sqrt(w)*27.211399)
```

