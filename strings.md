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

Josh has defined a slater determinant as 0b1101 to represent particles at 0&#593;1&#593;1&#946;, that is the lower states are to the right. We shall write this as '1011' that is lower states to the left. We are not constrained by the little-endian storage of binary data. We could store both the un-excited and excited determinants on the same string as say '1101:1011' but for now we'll save as two separate strings.

Let's take two determinants (Josh's single excitation example) 0b111 and 0b1000101. In our notation these will be '111' and '1010001'</br>
Firstly, lets get strings to same shape '111' -> '1110000'

```python
def normalise(det1, det2):
    #pad strings to equal length
    
    ld = max(det1, det2)
    if len(det1) != ld:
        det1 = det1 + '0' * (ld - len(det1))
    if len(det2) != ld:
        det2 = det2 + '0' * (ld - len(det2))
        
    return det1, det2   
```
Now we need to find out how many excitations our determinants represent. Let's put our determinants on top of each other</br>
'1110000'</br>
'1010001'</br>
Everywhere a '1' has become a zero represents an excitation (and where a '0' has become a '1'). The truth table for XoR satisfies this requirement for a bit implementation. We will define
```python
def excitations(det1, det2):
    #compute the total excitations
    
    excite = 0
    for i in range(len(det1)):
        if det1[i] != det2[i] : excite += 1
    
    return excite//2
```

Now we need to find the excitation levels. We can see (for single excitation) we need to find where the (first) '1' in determinant 1 has become a '0', then where the (first) '0' has become a '1' in determinant 2. For double excitations we just need to find the second occurrence of the preceding rule too.
```python
def levels(det1, det2):
    #compute the jumps
    
    ld = len(det1)
    jumps = []
    
    for i in range(ld):
        if det1[i] == '1' and det2[i] == '0':
            #hole has appeared
            for j in range(i+1, ld):
                if det2[j] == '1' and dt1[j] == '0':
                #particle has appeared               
                    jumps.append([i,j])
                    #set excited state to zero so don't count again
                    det2 = det2[:j] + '0' + det2[l+1:]
                    break
    
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
def phase(det1, det2):
    'compute the phase between the determinants
    
    ld = len(det1)
    
    for i in range(ld):
        if det1[i] == '1' and det2 == '0':
        #hole has appeared
            for j in range(i+1, ld):
                if det2[j] == '1' and dt1[j] == '0':
                #particle has appeared - get occupied between
                #hole and particle
                occupied = 0
                for k in range(i+1,j):
                    if det2[k] == '1': 
                        occupied += 1
                        #is this particle an excitation
                        if det1[k] == '0':
                            det2[k] = det2[:k] + '0' + det2[k+1:]
                            
      return -1**occupied
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
def combinations(m, k):
    #compute number of combinations of n things taken k at a time
    
    def fact(n):
        #compute factorial
        
        f == 1
        if n < 0:
            return -1
        elif n == 0:
            return 1
        else:
            for i in range(1,n + 1):
                f *= i
        return f
        
    return  fact(m)/(fact(k)*fact(m-k))  
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
first single 1&#946 excitations</br>
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
For water ((0, 1, 2, 3, 4, 5, 7, 10, 11, 20) state would be 111111010011000000001. Our code for this would be
```python
def binaryString(comb):
    #compute a binary string representation of combination as list
    
    sbinary = ''
    for i in range(max(comb)+1):
        if i in comb: sbinary = sbinary + '1'
        else: sbinary = sbinary + '0'
        
    return sbinary
```

Let's put together what we have...
```python
#molecule is Hydrogen in 3-21g basis
nElectrons = 2
nSpinOrbitals = nElectrons * 2

#how many excited determinants
determinants = combinations(nSpinOrbitals, nElectrons)

#get list of determinants
combinations = []
combinations = combinationList(combinations, '', 0, nSpinOrbitals-1, nElectrons-1)

#get binary string representation
binary = []
for det in range(combinations):
    binary.append(binaryString(combination[det]))
```
We now have to build our excited Hamiltonian. This will have dimension len(binary) x len(binary)</br>
The code will look like
```python
for i in range(len(binary)):
    for j in range(i+1):
        det1 = binary[i]
        det2 = binary[j]
        element = hamiltonianElement(det1, det2)
        H[i,j] = H[j,i] = element
```
Now to define hamiltonianElement. Firstly it must call *normalise* to make determinants equal length, then *excitations* to get the degree of the excitation, then if degree is <= 2 continue to process. It must now call *levels* to get the excitations themselves and finally call *phase* for the phase. </br>
If it's a double excitation we will have 4 values as say, \[0,3] and \[1,5] this is interpreted as the phase times <01||35>. </br>
For a single excitation say,\[1,6] this will be phase times &#931; <1n||6n>. Where n are common elements between the determinants. We need to calculate the common elements now. By elements we mean particle so its really just a logical AND. So
```python
def commonStates(det1, det2):
    #compute common states between determinants
    ld = min(len(det1), len(det2))
    
    common = []
    for i in range(ld):
        if det1[i] == '1' and det2[i] = '1': common.append(i)
        
    return common
```
there is also a one-body contribution which for \[1,6] would be H<sub>sc</sub>\[1,6], where H<sub>sc</sub> is the molecular spin core Hamltonian. Finally the zero degree exitations. These are m are the common states (namely anywhere there is a '1' in either determinant) so there is a H<sub>sc</sub>\[m,m] contribution and if m=n there is a pahse times a half times &#931; <mn||mn> for all combinations ie
```python
n = commonStates(det1, det2)

sum = 0.0
for i in range(n):
    for j in range(n):
        sum += <i,i||i,j>
        
sum += 0.5 * phase
```
That's how to build the excited Hamiltonian. To build the molecular basis spin core Hamiltonian firstly bring core Hamiltonian to molecular basis C<sup>T</sup>.H.C, then take kronecker product with I<sub>2</sub> ie H<sub>sc</sub> &#8855; I<sub>2</sub>, where I<sub>2</sub> is the 2x2 identity tensor.






    



  
 
