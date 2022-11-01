**Excited Determinants**

Josh pointed me in the direction of [this paper](https://hal.archives-ouvertes.fr/hal-01539072/document). It uses bit manipulation to implement Slater-Condon rules. The Slater-Condon rules are:
The essential results of the Slater-Condon rules are:

1.   The full N! terms that arise in the N-electron Slater determinants do not have to be treated explicitly, nor do the N!(N! + 1)/2 Hamiltonian matrix elements among
     the N! terms of one Slater determinant and the N! terms of the same or another Slater determinant
2.   All such matrix elements, for any one- and/or two-electron operator can be expressed in terms of one- or two-electron integrals over the spin-orbitals that appear 
     in the determinants.
3.   The integrals over orbitals are three or six dimensional integrals, regardless of how many electrons N there are.
4.   These integrals over mo's can, through the LCAO-MO expansion, ultimately be expressed in terms of one- and two-electron integrals over the primitive atomic 
     orbitals. It is only these ao-based integrals that can be evaluated explicitly (on high speed computers for all but the smallest systems).

For single excitations:</br>
  <&#916;|&#937;<sub>1</sub>|&#916;<sub>i</sub><sup>j</sup>> = <&#981;<sub>i</sub>|&#937;<sub>1</sub>|&#981;<sub>j</sub>> , where &#937;<sub>n</sub> is an n-body operator</br>
  <&#916;|&#937;<sub>2</sub>|&#916;<sub>i</sub><sup>j</sup>> = &#931; <&#981;<sub>i</sub>&#981;<sub>k</sub>|&#937;<sub>2</sub>|&#981;<sub>j</sub>&#981;<sub>k</sub>>  - <&#981;<sub>i</sub>&#981;<sub>k</sub>|&#937;<sub>2</sub>|&#981;<sub>k</sub>&#981;<sub>j</sub>> 
  
For double excitations:</br>
  <&#916;|&#937;<sub>1</sub>|&#916;<sub>ik</sub><sup>jl</sup>> = 0</br>
  <&#916;|&#937;<sub>2</sub>|&#916;<sub>ik</sub><sup>jl</sup>> = &#931; <&#981;<sub>i</sub>&#981;<sub>k</sub>|&#937;<sub>2</sub>|&#981;<sub>j</sub>&#981;<sub>l</sub>>  - <&#981;<sub>i</sub>&#981;<sub>k</sub>|&#937;<sub>2</sub>|&#981;<sub>l</sub>&#981;<sub>j</sub>> 
  
All higher excitations are 0. The paper states to implement the rules you need</br>
1.  to find the number of spin-orbital substitutions between two determinants
2.  to  find  which  spin-orbitals  are  involved  in  the substitution
3.  to  compute  the  phase  factor  if  a  reordering  of the spin-orbitals has occured.

The paper uses bit manipulation which is certainly efficient, but does it need to be? Surely the evaluation of the integrals will be the critical stage computationally. We could probably do the 3 steps above using strings. Python strings are limited in length only by the machine RAM so you do not need to worry about stitching variables together when 64 bits are used as in the bit method. Let's see how we go...

Josh has defined a slater determinant as 0b1101 to represent particles at 0&#593;1&#593;1&#946;, that is the lower states are to the right. We shall write this as '1011' that is lower states to the left as this seems more natural. We are not constrained by the little-endian storage of binary data. We could store both the un-excited and excited determinants on the same string as say '1101:1011' but for now we'll save as two separate strings.

Let's take two determinants (Josh's single excitation example) 0b111 and 0b1000101. In our notation these will be '111' and '1010001'</br>
Firstly, lets get strings to same shape '111' -> '1110000'

```python
def normalise(da, db):
    #pad strings to equal length
    
    ld = max(da, db)
    if len(det1) != ld:
        da = da + '0' * (ld - len(da))
    if len(db) != ld:
        db = db + '0' * (ld - len(db))
        
    return da, db  
```
Now we need to find out how many excitations our determinants represent. Let's put our determinants on top of each other</br>
'1110000'</br>
'1010001'</br>
Everywhere a '1' has become a zero represents an excitation (and where a '0' has become a '1').  We will define
```python
def excitations(da, db):
    #compute the total excitations
    
    excite = 0
    for i in range(len(da)):
        if (da[i] == '1') and (db[i] == '0' : excite += 1
    
    return excite
```

Now we need to find the excitation levels. We can see (for single excitation) we need to find where the (first) '1' in determinant 1 has become a '0', then where the (first) '0' has become a '1' in determinant 2. For double excitations we just need to find the second occurrence of the preceding rule too. Once we have found the excitation in the second determinant we must remove it to stop counting it twice, this means we must save a copy of second determinant to restore on exit.
```python
def levels(da, db):
    #compute the jumps
    
    ld = len(da)
    jumps = []
    t = db

    for i in range(ld):
        if da[i] == '1' and db[i] == '0':
            #hole has appeared
            for j in range(ld):
                if db[j] == '1' and da[j] == '0':
                #particle has appeared               
                    jumps.append([i,j])
                    #set excited state to zero so don't count again
                    db = db[:j] + '0' + db[j+1:]
                    break
    db = t

    return jumps
```
Now for the phase. The paper states *The phase is calculated as −1<sup>Nperm</sup>, where Nperm is the number permutations necessary to bring the spin-orbitals on which the holes are made to the positions of the particles. This number is equal to the number of occupied spin-orbitals between these two positions.*</br> So let's look at our example,</br>
'1110000'</br>
'1010001'</br>
The number of occupied spin orbitals between where the hole is (1) and where the particle is (6) is 1 (position 2) so the phase is -1<sup>1</sup> = -1. Let's look at the second example Josh gives ie 0b111 -> 0b101010, which we would write '111000' -> '010101'. This is obvoiusly a double excitation (0->3 and 2->5).</br> The phase is </br>
'111000'</br>
'010101'</br>
As above the first hole is (0) and the first particle is (3) and there is 1 occupied orbital between them (1), but there is a second hole at (2) and a second particle at (5) with 1 occupied orbital between them (3) - **but this is an excitation**. The paper notes *For a double excitation, if the realization of the first excitation introduces a new orbital between the hole and the particle of the second excitation  (crossing of the two excitations), an additional permutation is needed...* 
So we have 1+1+1=3 and the phase is -1<sup>3</sup> = -1. </br>
Since n+1 has the same parity n-1 (odd or even) instead of adding an extra permutation we could just not count the excited state in the first place. So we would say we have a hole at (0) and an excitation at (3) with one occupied state between, (3) now becomes '0' as we've dealt with it. We have a second hole at (2) and a particle at (5) with no occupied states between, so permutations are 1 and phase is -1.

```python
def phase(da, db):
    #compute the phase between the determinants
    
    t = db
    ld = len(da)
    occupied = 0

    for i in range(ld):
        if da[i] == '1' and db[i] == '0':
        #hole has appeared
            for j in range(ld):
                if db[j] == '1' and da[j] == '0':
                    #particle has appeared - get occupied between hole and particle
                    m = min(i,j)
                    n = max(j,i)
                    occupied += db.count('1', m+1, n)

                    db = db[:j] + '0' + db[j+1:]
                    break
     
    db = t  

    return pow(-1,occupied)
```
How can we get the permutations of the determinants - use scipy.special.comb to get the total number of combinations. As an example H<sub>2</sub> in 3-21g basis
there are 2 electrons in 4 basis functions, so
```python
print('For hydrogen molecule in 3-21g basis with 2 electrons and 8 (spin) basis functions there are.')
print(scipy.special.comb(2,8), ' determinants')
```
which gives an answer of 28 ie 8 things taken 2 at a time <sup>8</sup>C<sub>2</sub> = !8/!2 !(8-2) = 8.7/2 = 28
</br>Actually we could calculate our own combinations...
```python
def determinantCount(m, k):
    #compute number of combinations of n things taken k at a time
    
    def fact(n):
        #compute factorial
        
        f = 1
        if n < 0:
            return -1
        elif n == 0:
            return 1
        else:
            for i in range(1,n + 1):
                f *= i
        return f
        
    return  int(fact(m)/(fact(k)*fact(m-k)))
```
what are those combinations? </br>
(0,7)(0,6)(0,5)(0,4)(0,3)(0,2)(0,1)</br>
(1,7)(1,6)(1,5)(1,4)(1,3)(1,2)</br>
(2,7)(2,6)(2,5)(2,4)(2,3)</br>
(3,7)(3,6)(3,5)(3,4)</br>
(4,7)(4,6)(4,5)</br>
(5,7)(5,6)</br>
(6,7)</br>
this time itertools comes to the rescue
```python
for comblist in itertools.combinations(range(2), 8):
```
will generate the list in order (0,1)->(0,2)...->(6,7)

Here are the corresponding states...</br>
first single 1&#946; excitations</br>
(0,1) = '11'</br>
(0,2) = '101'</br>
(0,3) = '1001'</br>
(0,4) = '10001'</br>
(0,5) = '100001'</br>
(0,6) = '1000001'</br>
(0,7) = '10000001'</br>
then single 1&#945;</br>
(1,2) = '011'</br>
(1,3) = '0101'</br>
(1,4) = '01001'</br>
(1,5) = '010001'</br>
(1,6) = '0100001'</br>
(1,7) = '01000001'</br>
now multiple excitations</br>
(2,3) = '0011'</br>
(2,4) = '00101'</br>
(2,5) = '001001'</br>
(2,6) = '0010001'</br>
(2,7) = '00100001'</br>
(3,4) = '00011'</br>
(3,5) = '000101'</br>
(3,6) = '0001001'</br>
(3,7) = '00010001'</br>
(4,5) = '000011'</br>
(4,6) = '0000101'</br>
(4,7) = '00001001'</br>
(5,6) = '0000011'</br>
(5,7) = '00000101'</br>
(6,7) = '00000011'</br>

We could calculate our own list...</br>
```python
def combinationList(combs, group, start, stop, level):
    #compute the combinations for taking n things k at a time
    
    for i in range(start, stop+1):

        if level == 0: 
            s = (group + ',' + str(i))[1:]
            combs.append(list(map(int, s.split(','))))
        
        combinationList(combs, group + ',' + str(i), i+1, stop, level-1)
    
    return combs

#for our problem run as...
combs = []
n = 8
k = 2
combinations = combinationList(combs, '', 0, n-1, k-1)
```

How do we generate the binary strings? </br>
Well (0,1) we say is 2<sup>0</sup>+2<sup>1</sup> = 3 -> '11', as (1,3) would be 2<sup>1</sup>+2<sup>3</sup> = 2+8 = 10 = '0101' ie a '1' in positions 1 and 3</br>
For water ((0, 1, 2, 3, 4, 5, 7, 10, 11, 13) state would be 11111101001101. Our code for this would be
```python
def binaryString(comb, nOrbitals):
    #compute a binary string representation of combination as list
    
    sbinary = ''
    for i in range(max(comb)+1):
        if i in comb: sbinary = sbinary + '1'
        else: sbinary = sbinary + '0'
        
    return sbinary + '0' * (nOrbitals - len(sbinary))
```

Let's put together what we have...
```python
#molecule is Hydrogen in 3-21g basis
nElectrons = 2
nSpinOrbitals = nElectrons * 2

#how many excited determinants
determinants = combinationCount(nSpinOrbitals, nElectrons)

#get list of determinants
combinations = []
combinations = combinationList(combinations, '', 0, nSpinOrbitals-1, nElectrons-1)

#get binary string representation
binary = []
for det in range(combinations):
    binary.append(binaryString(combination[det]))
binary.sort()
```
We now have to build our excited Hamiltonian. This will have dimension len(binary) x len(binary)</br>
The code will look like
```python
def buildFCIhamiltonian(determinants, eriMOspin, Hp):
    #compute the full FCI Hamiltonian

    nH = len(determinants)
    fciH = np.zeros((nH, nH))

    for i in range(len(determinants)):
        for j in range(0, i+1):

            da = determinants[i]
            db = determinants[j]

            element = hamiltonianElement(da, db, eriMOspin, Hp)
            fciH[i,j] = fciH[j,i] = element

    return fciH
```
Now to define hamiltonianElement. Firstly it must call *excitations* to get the degree of the excitation, then if degree is <= 2 continue to process. It must now call *levels* to get the excitations themselves and finally call *phase* for the phase. </br>
If it's a double excitation we will have 4 values as say, \[0,3] and \[1,5] this is interpreted as the phase times <01||35>. </br>
For a single excitation say,\[1,6] this will be phase times &#931; <1n||6n>. Where n are common elements between the determinants. We need to calculate the common elements now. By elements we mean particle so its really just a logical AND. So
```python
def commonStates(da, db):
    #compute common states between determinants
    ld = min(len(da), len(db))
    
    common = []
    for i in range(ld):
        if da[i] == '1' and db[i] == '1': common.append(i)
        
    return common
```
For single excitations there is also a one-body contribution which for \[1,6] would be H<sub>sc</sub>\[1,6], where H<sub>sc</sub> is the molecular spin core Hamltonian. Finally the zero degree exitations. These are m are the common states (namely anywhere there is a '1' in either determinant) so there is a H<sub>sc</sub>\[m,m] contribution and if m=n there is a phase times a half times &#931; <mn||mn> for all combinations giving overall
```python
def hamiltonianElement(da, db, eriMOspin, Hp):
    #compute an individual Hamiltonian element

    #get number of excitions
    excites = excitations(da, db)

    #excitions above 2 do not contribute
    if excites > 2 : return 0.0

    #get phase and excitations for 2,1 and 0
    theta = phase(da,db)
    jump = levels(da, db)
    if excites == 2:
        return theta * eriMOspin[jump[0][0], jump[1][0], jump[0][1], jump[1][1]]

    elif excites ==1:
        common = commonStates(da, db)
        f = Hp[jump[0][0],jump[0][1]]
        for i in common:
            f += eriMOspin[jump[0][0], i, jump[0][1], i]

        return theta * f

    elif excites == 0:
        #really only 1 determinant just as cheap to use common
        common = commonStates(da, db)
        f = 0.0
        for i in common:
            f += Hp[i, i]
        for i in common:
            for j in common:
             f += 0.5 * eriMOspin[i, j, i, j]

        return theta * f
```
That's how to build the excited Hamiltonian. To build the molecular basis spin core Hamiltonian firstly bring core Hamiltonian to molecular basis C<sup>T</sup>.H.C, then take kronecker product with I<sub>2</sub> ie H<sub>sc</sub> &#8855; I<sub>2</sub>, where I<sub>2</sub> is the 2x2 identity tensor. Alternatively, if you don't have a Kron function you could do...
```python
def buildCoreMOspin(spinOrbitals, eigenVectors, core):
	#transform core -> MO -> spin basis

	corespin = zeros((spinOrbitals, spinOrbitals))
	eigenspin = zeros((spinOrbitals, spinOrbitals))

	for p in range(0, spinOrbitals):
		for q in range(0, spinOrbitals):
			corespin[p,q]  = fock[int(p/2), int(q/2)] * ((p % 2) == (q % 2))
			eigenspin[p,q] = eigenVectors[int(p/2), int(q/2)] * ((p % 2) == (q % 2))

	return dot(eigenspin.T, dot(corespin, eigenspin))
```

For CISD we need a different approach. For H<sub>2</sub> we have two parts, the first part is the nElectron ground state eg '11'. Aside from the ground state we want all single '0' substitutions ie '01' and '10', then all double substitutions ie '00'. The second part is the nOrbitals - nElectron part. This must be all combinations of the nOrbitals - nElectron taken k at a time, where k is a number to make total electron count nElectrons. Let's look at an example, \
For H<sub>2</sub>, \
'11' + '000000' \
'01', '10' + '000001', '000010', '000100', '001000', '010000', '100000' \
'00' + '000011', '000101', '001001', '010001', '100001', '000110', '001010', '010010', '100010', '001100', '010100', '100100', '011000', '101000' \
total 28. We have code to do this
```python
def configurations(nElectrons, nOrbitals, type = 'S'):

    determinants = []
    pad = nOrbitals - nElectrons

    def subDeterminant(n, k, bit):
        #components of full determinant

        sub = []
        comb = []
        combinationList(comb, '', 0, n-1, k)

        for i in comb:
            s = bit[0] * n
            for j in range(k+1):
                s = s[:i[j]] + bit[1] + s[i[j]+1:]
            sub.append(s)

        return sub

    #groundstate
    if 'G' in type:
        determinants.append('1' * nElectrons + '0' * pad)

    #generate groundstate single excitations
    if ('S' in type) and (nElectrons > 0):

        pre = subDeterminant(nElectrons, 0, '10')
        post = subDeterminant(pad, 0, '01')

        for i in pre:
            for j in post:
                determinants.append(i+j)
        
    #generate groundstate double excitations
    if ('D' in type) and (nElectrons > 1):

        pre = subDeterminant(nElectrons, 1, '10')
        post = subDeterminant(pad, 1, '01')

        for i in pre:
            for j in post:
                determinants.append(i+j)

    #generate groundstate triple excitations
    if ('T' in type) and (nElectrons > 2):

        pre = subDeterminant(nElectrons, 2, '10')
        post = subDeterminant(pad, 2, '01')

        for i in pre:
            for j in post:
                determinants.append(i+j)

    #generate groundstate quadruples excitations
    if ('Q' in type) and (nElectrons > 3):

        pre = subDeterminant(nElectrons, 3, '10')
        post = subDeterminant(pad, 3, '01')

        for i in pre:
            for j in post:
                determinants.append(i+j)

    determinants.sort()

    return determinants
```
We've implemented configurations( 'S' - for CIS and 'GSD' for CISD.

If we want to get back the alpha and beta spin determinants we can use
```python
def getSpinState(da, nBasis):
    #get the alpha and beta spin states

    alpha = da[0::2]
    beta  = da[1::2]

    return  alpha + '0'*(nBasis - len(alpha)) , beta + '0'*(nBasis - len(beta))
    
print(getSpinState('00100101001111',7))

('0100011', '0011011')
```
The occupancy of the spatial can be found by
```python
def getOccupancy(determinant, nBasis):
    #get the unoccupied, singe and double occupancy (spatial) orbitals

    alpha , beta = getSpinState(determinant, nBasis)

    occupancy = ''
    for i in range(len(alpha)):
        occupancy += str(int(alpha[i]) + int(beta[i]))

    return occupancy

print(getOccupancy('00100101001111',7))

0111022

```
### Bit Manipulation
Using textual bit strings makes thing clear, and for our needs are fast enough plus they can work with any number of orbitals. 
However, it has been pointed out to me that in eg selected Configuration Interaction determinants need to be generated in a 'as needed'
basis, rather than in one operation as is done here. So we will look at the more traditional bit-opertions techniques limiting to 64bits (spin orbitals).
A python integer word is 64 bits long, so does that mean we are limited to bit strings of 64 orbitals...
```python
h = '0b1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111'
print('Decimal representation...',int(h,2))

1Decimal representation...267650600228229401496703205375
```
this is a bit string of 100 bits and python handles it OK. Josh Goings points out that with numpy use dtype=object 
to prevent it using int64 and overflow occuring

The number of combinations can be found by using scipy.special.comb, so nOrbitals filled by nElectrons at a time is (eg water)
```python
nBasis = 7
spinOrbitals = nBasis * 2
nElectrons = 10

from scipy.special import comb
print("Number determinants... ",comb(spinOrbitals, nElectrons, exact=True))

Number of determinants...1001
```
Now to get the combinations of where the nElectrons are in the spinOrbitals we can use
```python
from itertools import combinations
for determinant in combinations(range(spinOrbitals), nElectrons):
    print(determinant)
```
this returns values in determinant as eg (3, 4, 5, 6, 7, 8, 9, 11, 12, 13), that is nElectrons positions in a 0-spinOrbitals-1 range.
We now need this as a bit string which for the above example would be 00011111110111. We are in a little-endian world though so we need to reverse
string to 11101111111000 and add '0b' in front, then we will encode it bit string as a denary digit
```python
determinant = (3, 4, 5, 6, 7, 8, 9, 11, 12, 13)

def bString(determinant, spinOrbitals):
    #convert determinat to a bit string integer

    bits = '0b'
    for i in range(spinOrbitals-1, -1, -1): 
        if i in determinant: bits += '1'
        else: bits += '0'

    return int(bits,2)

print('Bit string representation...',bString(determinant, spinOrbitals))

Bit string representation...15352

```
Now put it all together to get a list of all combinations
```python
def bCombinationList(spinOrbitals, nElectrons):
    #create list of all combinations of nElectrons in spinOrbital orbitals

    determinantList = []

    for determinant in combinations(range(spinOrbitals), nElectrons):

        determinantList.append(bString(determinant, spinOrbitals))

    return determinantList
```
we can get the holes (zeros) in excitaion da->db, but first we need to define two service routines
```python
#first non-zero trailing bit
def bRightZeros(n):
    #compute the number of rightmost zero bits

    return (n & -n).bit_length() - 1

#clear nth bit
def bSetZero(da, n):
    #set the nth bit to zero 0-based

    return ~(1 << n) & da

h = '0b1111111111111'
i = bSetZero(int(h,2), 4)
print('Bit 4 set to 0', bin(i))

Bit 4 set to 0 0b1111111101111
```
we can now find all the holes and particles
```python
def occupancy(da, db, type = 'h'):
    #compute the holes between determinants

    k = 0
    hole = []
    if type == 'h': h = (da ^ db) & da
    else:  h = (da ^ db) & db

    while h != 0:
        p = bRightZeros(h)
        hole.append(p)

        h = bSetZero(h, p)
        k += 1

    return hole

da = 7
db = 42
h = occupancy(da, db, 'h')
print('Holes in 111->101010 ', h)

Holes in 111->101010  [0, 2]

p = occupancy(da, db, 'p')
print('Particles in 111->101010 ', p)

Particles in 111->101010  [3, 5]
```
we now need bitstring versions of our Hamiltonian builds. Actually the routine buildFCIhamiltonian can be used unchanged as it is just a 
conduit for the determinant list though to the hamiltonianElement routine. First the excitations, first how many are there? The degree of excitation is equivalent to the number of holes (zeros) created in first determinant of a pair or to the number of particles (ones) created in the second. The following routine will do this by xor'ing the two determinants and counting the '1''s ie where a 01 or 10 has occurred. Finally we must divide by 2 because we have counted both 01 and 10.
```python
def bExcitations(da, db):
    #the number of excitation between the two determinants

    return (bin(da ^ db).count('1')) >> 1

#test for '0b111' (7) and '0b101010' (42)
da = 7
db = 42
excitation = bExcitations(7, 42)
print('Number of excitations...', excitation)

Number of excitations... 2
```
which is correct as we have 2->5, 0->3. Now we need to know what the excitations are, but only if there are there are two or fewer. For **single** excitations...\
if da and db are the same there is no excitation \
it is more efficient to recode the holes and particles rather than using the routines we already wrote
```python
ab = da ^ db
hole     = ab & da
particle = ab & db
```
we define an array excite\[2,2] where \[i,.] for i=0 is degree of excitation and i=1 is position of orbital, \[.,i] for i=0 is hole and i=1 is particle.
```python
excite = np.zeros((2,2), dtype=object)

excite[0,0] = 1
excite[1,0] = bRightZeros(hole)
excite[0,1] = 1
excite[1,1] = bRightZeros(particle)
 ```
Now we need the phase defined as above
```python
lo = min(excite[1,0], excite[1,1])
hi  = max(excite[1,0], excite[1,1])

nPerm = bin(da &  (~(1 << lo+1)+1  & (1 << hi)-1)).count('1')
```
How does the last statement work? (1<<low+1)  makes 000...1**0**00...000 where the hole position is bold, ~ make this 111...0**1**11...111, and +1 makes 111...1**0**00...000 - so this is a mask with all bits greater than the hole are 1 and all bits at and less than the hole are 0. (1<<hi) sets the particle bit to 1 and '-1' sets a mask 111.... where the mask starts 1 bit to the right of the particle bit and'ing the two masks creates a mask 00...01111...1000...000, so all the bits between hole and particle are 1 and all others are 0. Now 'and' with the da and count the 1's for the number of 1's between the hole and particle. \
The phase is −1<sup>Nperm</sup>. If we 'and' an odd number with 1 we will get 1 or if we 'and' an even number with 1 we will get 0 so the phase can be written as
```python
phases = [1, -1]
phase = phases[nPerm & 1]
```
Now for the double excitations. Again we get the hole and particle->get their position->increase the number of excitations and store->store the position->clear the current bit->loop back to find a new position until none found.
```python
excite = np.zeros((3,2), dtype=object)

ab = da ^ db
hole     = ab & da
particle = ab & db

#the holes
iHole = 0
while hole != 0:
    iHole += 1
    excite[iHole,0] = bRightZeros(hole)
    excite[0,0] += 1
    
    hole = hole & (hole - 1)

#the particles
iParticle = 0
while particle != 0:
    iParticle += 1
    excite[iParticle,1] = bRightZeros(particle)
    excite[0,1] += 1
    
    particle = particle & (particle - 1)
    
#phase
for excitation in [1,2]:
    lo = min(excite[excitation,0],excite[excitation,1])
    hi = max(excite[excitation,0],excite[excitation,1])

    nPerm += bin(da &  (~(1 << lo+1)+1  & (1 << hi)-1)).count('1')

#orbital crossings
i = min(excite[1,0], excite[1,1])
a = max(excite[1,0], excite[1,1])
j = min(excite[2,0], excite[2,1])
b = max(excite[2,0], excite[2,1])
if j>i and j<a and b>a : nPerm += 1

phases = [1,-1]
phase = phases[nPerm & 1]
```
Next we will need a routine to find the common indices. Here we 'and' the two determinants which gives us a 1 where there are common elements and 0 elsewhere.\
We then create masks with a 1 in every possible position and if the 'and' with that mask has a 1 then that is a common element position.
```python
common = []
ab = da & db

for i in range(aandb.bit_length()):
    idx = ab & (1 << i)
    if idx: common.append(bRightZeros(1<<i))
```
This is all we need to do an fci computation using bits.\
Residues.  The residues generated from a reference determinant are a set of determinants that is generated by removing two electrons in all possible ways from the reference determinant. If we make two binary expressions that have just 1 bit each set to 1 in different positions, form to xor and then negate we will have a mask of all 1's except for two bits which are cleared. Now 'and' this with the original determinant (r) and the number of 1's are are the the number of occupied orbitals other than the two original bits in the xor mask (n). If the total occupancy of the original determinant minus n is 2 then particles in the two original bits can excite and r is a residue. Do this for all pairs of bits. Residues are important in selected CI methods as the residue array method.
```python
residues = []
occupancy = bin(da).count('1')

for i in range(spinOrbitals):
    mask1 = 1 << i
    for j in range(i):
        mask2 = 1 << j
	reducedOccupancy = da & ~(mask1 ^ mask2)
	
	if bin(reducedOccupancy).count('1') == (occupancy - 2):
	    residues.append(reducedOccupancy)
```	    
Next we have to find all ways of filling the holes in a residue. The technique is similar to the one for finding the residues. We form two one bit masks again, test if the position of a mask is 0 so a particle can move there. If determinant 'and' mask is 0 then there is hole and we can find this, because 0 = False, so not bool(da & mask) is the test. If we 'or' the mask with determinant we get the determinant with the hole set.
```python
determinants = []

for residue in residues:
    
    for i in range(spinOrbitals):
        mask1 = 1 << i
	if not bool(da & mask1):
	    p = da | mask1
	    
	    for j in range(i):
	        mask2 = 1 << j
		if not bool(p & mask2):
		    q = p | mask2
		    determinants.append(q)
```
The results are identical to the text bit strings. 
		    

    



  
 
