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

### alpha and beta eigenvectors
To get these we need to do a TDHF calculation. To avoid the left and right eigenvector problem for the moment we'll use the Hermitian form of the equations \
The basic form of the equations \
**A** = &#948;<sub>𝑖𝑗</sub>**𝑓**<sub>𝑎𝑏</sub> − &#948;<sub>𝑎𝑏</sub>**𝑓**<sub>𝑖𝑗</sub> + <bi||𝑗a> \
**B** = <bi||ja> \

(**A** - **B**)(**A** + **B**) |**X**<sub>n</sub> + **Y**<sub>n</sub>> =  	&#969;<sup>2</sup> |**X**<sub>n</sub> + **Y**<sub>n</sub>> (i) \
left multiply by (**A** - **B**)<sup>-&#189;</sup> and right multiply by (**A** - **B**)<sup>&#189;</sup> gives \
(**A** - **B**)<sup>&#189;</sup>**.**(**A** + **B**)**.**(**A** - **B**)<sup>&#189;</sup> **T**<sub>n</sub> = &#969;<sup>2</sup> **T**<sub>n</sub> (ii) \
writing (ii) as (left multiply by (**A** - **B**)<sup>&#189;</sup>) gives \
(**A** - **B**)(**A** + **B**) ((**A** - **B**)<sup>&#189;</sup>**T**<sub>n</sub>) =  &#969;<sup>2</sup> ((**A** - **B**)<sup>&#189;</sup> **T**<sub>n</sub>) \
comparing this Hermitian form with (i) gives  |**X**<sub>n</sub> + **Y**<sub>n</sub>>  =  ((**A** - **B**)<sup>&#189;</sup> **T**<sub>n</sub>) (iii) \
There is a corresponding Hermitian form developed from \
(**A** + **B**)(**A** - **B**) |**X**<sub>n</sub> - **Y**<sub>n</sub>> =  	&#969;<sup>2</sup> |**X**<sub>n</sub> - **Y**<sub>n</sub>> (iv) \
(**A** + **B**)<sup>&#189;</sup>**.**(**A** - **B**)**.**(**A** + **B**)<sup>&#189;</sup> **S**<sub>n</sub> = &#969;<sup>2</sup> **S**<sub>n</sub> (v) \
which gives |**X**<sub>n</sub> - **Y**<sub>n</sub>>  =  ((**A** + **B**)<sup>&#189;</sup> **S**<sub>n</sub>) (vi) \
Again |**X**<sub>n</sub> - **Y**<sub>n</sub>> and |**X**<sub>n</sub> + **Y**<sub>n</sub>> should be biorthonormal.

It is not necessary to separately solve for |**X**<sub>n</sub> - **Y**<sub>n</sub>> if you have |**X**<sub>n</sub> + **Y**<sub>n</sub>> as they are connected through \
(**A**+**B**)|**X**<sub>n</sub> + **Y**<sub>n</sub>> = &#969; |**X**<sub>n</sub> - **Y**<sub>n</sub>> or \
(**A**-**B**)|**X**<sub>n</sub> - **Y**<sub>n</sub>> = &#969; |**X**<sub>n</sub> + **Y**<sub>n</sub>> 


Running a TDHF calculation using the Hermitian form gives eigenvalues 
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

The eigenvectors from psi4 for alpha eigenvectors right for the first and second singles singlets - alpha and beta are the same. (It is clear the eigenvalues are not themselves normalised.)
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
The harpy implementation is 
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
We also have a gram-schmidt routine
```python
def gs(X, row_vecs=True, norm = False):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T
```
Need to transform the excitation eigenvalue into a No x Nv block of alpha (or beta) values
```python
alpha = []
for i in range(0,2*No, Nv):
    j = 2*Nv*i
    alpha.append(t[j:j+2*Nv:2])
alpha = np.array(alpha).reshape(No,Nv)
```
This gives us
```
[[-3.73839867e-19  7.72362896e-17]
 [ 3.22297056e-17  3.88439857e-15]
 [ 5.24681720e-14  3.17527244e-17]
 [-2.06926883e-16 -7.15164108e-14]
 [-4.01159972e-01  6.24817801e-14]]
 
[[-4.70155477e-20 -1.57571902e-17]
 [ 4.17835509e-17  8.42240232e-16]
 [-7.96646511e-14 -2.30714974e-16]
 [ 1.39987328e-16 -6.09804647e-14]
 [ 5.20446267e-14 -4.42170995e-01]]
```
These are not right.
### Calculation of properties
For the time being I'll work with the psi4 values to verify the procedure for calculating properties \
To calculate the transition dipole we evaluate (2)<sup>1/2</sup>**&#956;.**|**X**<sub>n</sub>+**Y**<sub>n</sub>>
```python
psi = [[ 1.24887231e-16,  5.76788541e-19],
       [-1.08892241e-15, -6.05807431e-17],
       [ 2.13542747e-15, -5.54872337e-16],
       [ 3.94632694e-15,  3.10090564e-16],
       [ 9.52476119e-01,  1.79097524e-15]]

tdm = []
for i in range(3):
    tdm.append(np.sqrt(2)* np.sum(psi*muMO[i]))
print('Transition Dipole ', tdm)

Transition Dipole [-4.1925798959575217e-16, 3.929517040051327e-15, 0.09454081962660157]
```
psi4 values are \[ 9.04458476e-16 -1.52641890e-16  9.45406317e-02]

Next we can get the oscillator strength
```python
os = (2/3) * np.sqrt(w[9]) * np.sum(np.array(tdm)**2)
print('Oscillator strength ', os)

Oscillator strength  0.0021139975452614036

```
psi4 value is OSCILLATOR STRENGTH (LEN) = 0.0021139901732160796

Going back to the equation (i) \
(**A** - **B**)(**A** + **B**) |**X**<sub>n</sub> + **Y**<sub>n</sub>> =  	&#969;<sup>2</sup> |**X**<sub>n</sub> + **Y**<sub>n</sub>> (i)\
The matrix (**A** - **B**)(**A** + **B**) is non-Hermitian, giving rise to left and right eigenvectors the right eigenvectors as above are |**X**<sub>n</sub> + **Y**<sub>n</sub>> and the left-eigenvectors are <**X**<sub>n</sub> - **Y**<sub>n</sub>|. \
We can get in left-eigenvectors from the right ones by the relation [given here](http://link.aip.org/link/doi/10.1063/1.477483?ver=pdfcov), which is \
(**A** + **B**) |**X**<sub>n</sub> + **Y**<sub>n</sub>> = &#969; |**X**<sub>n</sub> - **Y**<sub>n</sub>> \
This allows us to calculate velocity gauge properties. We can calculate left eigenvectors eg
```python
leftEigenvector = np.dot(A+B, rightEigenvector).dot(np.diag(1.0, np.sqrt(w)))
```
For reference the left alpha eigenvectors are
```
[[0.00000000e+00 0.00000000e+00]
 [0.00000000e+00 0.00000000e+00]
 [0.00000000e+00 0.00000000e+00]
 [6.09409487e-15 1.00710682e-15]
 [1.04989509e+00 1.64909657e-14]]
 
[[ 0.00000000e+00  0.00000000e+00]
 [-9.99179002e-16  6.62025083e-16]
 [-1.70694081e-16  1.89986020e-16]
 [ 4.46370254e-15 -2.81796845e-15]
 [-3.62245267e-14  1.03058747e+00]]
```
We could instead of <&#968;<sub>0</sub>|**&#956;**|&#968;<sub>i</sub> ie length gauge use <&#968;<sub>0</sub>|**&nabla;**|&#968;<sub>i</sub>> ie velocity gauge \
and calculate as (2)<sup>1/2</sup> *|**X**<sub>n</sub>-**Y**<sub>n</sub>>
ELECTRIC DIPOLE TRANSITION MOMENT (VEL) = \[7.01779000e-18 2.33393269e-16 1.37034011e-01]
OSCILLATOR STRENGTH (VEL) = 0.035286473435577954

The magnetic transition dipole in length gauge is calculated as (1/2)<sup>1/2</sup>**&omega;.**|**X**<sub>n</sub>-**Y**<sub>n</sub>> ,where **&omega;** are the angular momentum integrals.

The length gauge rotatory strength is calculated as 
```python
lgrs = np.sum(etdm_l * mtdm)_l
```
where etdm_l is the electronic transition dipole momemnt in the length gauge and mtdm_l is the magnetic transition dipole moment in the length gauge. \
The velocity gauge rotatory strength is calculated as 
```python
mgrs = -np.sum(etdm_v * mtdm_l)/w
```
For both one-photon absorption (OPA) and electronic circular dichroism (ECD) the poles are the excitation energies. The residues are: \
    The square of the transition electric dipole moment for OPA: |**&mu;** |<sup>2</sup>  
    The rotatory strength **R**<sub>ij</sub> for ECD. This is the imaginary part of the product of transition electric and magnetic dipole moments: Im(**&mu;**<sub>ij</sub>**.m**<sub>ij</sub>) 
The formulas are ![](formula.png)

**note** \
The Hermitian treatment means you can use symmetric solvers for the eigen problem. However, for the |**X**+**Y**> equation (iii) there is another Hermitian one for |**X**-**Y**> and once you've solved for these quantities you still need to biorthonormalise them to be correct eigenvectors to the problem? I now think that there is no analytic way of generating a biorthonormal set of eigenvectors from a general pair of eigenvectors and it must be done by an iterative process.

## Tamm-Dancoff Approximation
Set **B** = 0. For water eigenvalues
```
[ 0.28725552  0.28725552  0.28725552  0.34442501  0.34442501  0.34442501
  0.35646178  0.36598901  0.36598901  0.36598901  0.39451381  0.39451381
  0.39451381  0.41607175  0.5056283   0.51429     0.51429     0.51429
```
From psi4 using
```
res = tdscf_excitations(wfn, states=5, tda=True)

---> 0.3564619474642684
[[ 3.62517239e-17 -7.06179792e-20]
 [-3.74940352e-16 -5.73760059e-18]
 [ 6.57671825e-16 -9.83232306e-17]
 [ 2.26894887e-15  1.76536859e-17]
 [ 1.00000000e+00  1.23948800e-14]]
---> 0.4160718299870798
[[ 9.79766595e-17  2.21455727e-17]
 [ 2.23711109e-15 -3.39925229e-16]
 [-3.45280790e-16  6.72851332e-15]
 [ 1.26680223e-14  5.58514303e-15]
 [-1.23948800e-14  1.00000000e+00]]
---> 0.505628345495891
[[-7.84171394e-04 -2.34871529e-16]
 [-6.40451105e-02  2.66235698e-15]
 [-5.17607623e-14 -3.44767100e-01]
 [-9.36500537e-01 -7.29393875e-14]
 [ 2.06698855e-15  1.44451766e-14]]
---> 0.5551918484163628
[[-3.58582709e-17  1.08027457e-03]
 [-1.17579349e-14  1.50332457e-03]
 [ 5.52669710e-01 -8.26655554e-14]
 [-6.41864567e-14  8.33398323e-01]
 [-3.78179147e-16 -4.59046921e-15]]
```
From harpy for 3rd eigenvalues and after normalising the eigenvector I get
```
---> 0.5056283
[[ 7.84170883e-04  1.10912811e-16]
 [ 6.40451185e-02 -1.25496455e-15]
 [ 9.75763912e-15  3.44766688e-01]
 [ 9.36500688e-01  1.67716412e-14]
 [ 3.81388604e-30 -3.63440589e-16]]
```
upto a sign these agree. (dipole should have minus but still odd wrong sign)

The excitations, transition dipole (length) and oscillator strengths from psi4 for first 5 singlet states are
```
0.35646194746426896
[-1.00985303e-15  7.04945335e-16  9.92577450e-02]
0.0023412658233158954
0.41607182998707265
[-8.04597693e-15  7.82930516e-16  7.43339206e-16]
1.8280331676044886e-29
0.5056283454958909
[-7.77099088e-15  4.38874739e-01 -3.07772645e-16]
0.06492639995293956
0.5551918484163678
[-2.04425928e-01 -8.79344690e-14 -1.08171461e-16]
0.015467630166368004
0.6553184873358628
[1.69281940e+00 1.28985824e-14 1.58400793e-17]
1.2519368251364011
```
velocity gauge works but again a sign difference here and there (dipole integral signs agree so difference in the eigenvectors). Magnetic dipole next then find example with non-zero rotary strengths.
I can now reproduce these results (modulo a sign here and there). In TDA the left eigenvectors are equal to the right eigenvectors. I get
```
0.3564617754312789
[4.12046668e-15, 2.89530978e-21, -0.09925794]
0.0023412739974649588,
0.41607174996907337
[-1.10887617e-14, 3.92186240e-15, 1.95094676e-16]
3.83839238e-29
0.505628299588154
[1.20886205e-14, 0.43887432, 4.15792989e-31]
0.06492626973606519
0.5551918931220665
[-0.20442409, -3.19631675e-15, -3.27682611e-17]
0.015467354078846973
0.6553184580683947
[1.69281951, 1.14084334e-14, -3.29033947e-17]
1.2519369374141582
```










