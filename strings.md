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

    #get pahse and excitations for 2,1 and 0
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

For CISD we need a different approach, first get the residues - which are all combinations of removing two particles from a ground state determinant (nElectrons). 
```python
combs = []
nElectrons = 10
nOrder = 2
combinations = combinationList(combs, '', 0, nElectrons-1, nOrder-1)

def residues(combinations, n):
    #compute residues

    residueList = []

    for i in range(len(combinations)):
        residue = ''
        for j in range(n):
            s = '1'
            if j in combinations[i]: s = '0'
            residue = residue + s
        
        residueList.append(residue)

    return residueList
    
holes = residues(combinations, nElectrons)
```
Our *residues* strings have (nElectron - 2) particles. We now need to add 2 particles in (nOrbitals-nElectron) spaces. Eg </br>
For H<sub>2</sub> combinations will be (0,1) which will result in residueList being \['00']. For water our first residue will be '00111111'. First right pad to number of orbitals (spin) ie 14 '001111110000'. We will now have 3 cases. </br>
Case 1. Both zeros in unpadded residue are now '1's - this gives 1 case '11111111110000' </br>
Case 2. One zero in unpadded residue is a '1' - this, for water, gives right hand strings of '1000', '0100', '0010' and '0001' = 4 times 2 for other zero in unpadded residue. So total of 8. <sup>4</sup>C<<sub>3</sub> </br>
Case 3. No zeros in unpadded residue are '1's - this, for water, gives right hand strings of '1100', '1010', '1001', '0110', '0101' and '0011' = 6,So total <sup>4</sup>C<sub>2</sub> </br>
So a total of 15 for each residue. In general, there will be 1 + <sup>n</sup>C<sub>3</sub> + <sup>n</sup>C<sub>2</sub> ,where n is number of spin orbitals - the number of electrons. 

```python
def particles(residues, nElectrons, nOrbitals):
   #add particles back to spin orbital

    determinants = []
    n = nOrbitals - nElectrons

    #loop over the residues
    for i in range(len(residues)):

        #'11' in holes
        determinant = residues[i].replace('0', '1') + '0' * n
        determinants.append(determinant)

    
        #'01' or '10' in holes
        pad = []
        k = 1
        pad = combinationList(pad, '', 0, n-1, k-1)

        residuea = residues[i].replace('0', '1', 1)
        k = residues[i].rfind('0')
        residueb = residues[i][:k] + '1' + residues[i][k+1:]

        for j in range(len(pad)):
            s = binaryString(pad[j], n)
            determinants.append(residuea + s)
            determinants.append(residueb + s)
        
        #'00' in holes
        pad = []
        k = 2
        pad = combinationList(pad, '', 0, n-1, k-1)
        for j in range(len(pad)):
            s = binaryString(pad[j], n)
            determinants.append(residues[i] + s)
        
    unique = []
    for i in determinants:
        if i not in unique: unique.append(i)

    unique.sort()   
    return unique
```
Can we extend to do T and Q? Residues will produce determinants with as many holes as specified by nOrder in the combination passed to it. But particles is hardcoded for 2 particles. Say we have 4 electrons and 8 spin orbitals and we want triples so nOrder is 3. Calling combinationList(combs, '', 0, 3, 2) will produce combinations \[0, 1, 2], \[0, 1, 3], \[0, 2, 3], \[1, 2, 3]. We now need particles to produce from these residues - > '0001', '0010', '0100', '1000'. Now for say '0001', this will go to '00010000' (for 8 spin orbitals). We have nOrder cases \
0. Replace nOrder holes in original residue with '1's and do nothing with padded (nOrbitals - nElectrons)=k '0'. \
1. Replace nOrder-1 holes in original residue with '1's and replace k-1 '0's in padded with '1's. \
... \
n. replace nOrder-n holes in original residue with '1's and replace k-n '0's in padded with '1's. \
... \
nOrder.  replace nOrder-nOrder holes in original residue with '1's and replace k-nOrder '0's in padded with '1's. \

In our example, \
[0] '11110000' \
[1] '0111' + '1000','0100','0010', '0001' and '1011' + '1000','0100','0010', '0001' and '0111' + '1000','0100','0010', '0001'. \
[2] '0011' + '0011', '0101','0110', '1001', '1010', '1100' and '0101' + '0011', '0101','0110', '1001', '1010', '1100' and '1001' + '0011', '0101','0110', '1001', '1010', '1100' \
[3] '0000' + '1111'

In general, \









    



  
 
