## Cluster Operator Generation Using Sympy

Sympy has a module called secondquant and I want to use this to automatically generate some cluster operators.
+ Generate some symbols\
  Quantum chemistry is pretty heavy on labels, we have lots of tensors, so how using sympy do we define some labels to use with our tensors. Unlike in Python where we just just use the variables we want without declaration, in sympy we must declare our variables/labels
  
```python
  i, j = symbols('i, j', cls=Dummy)
  a, b, c = symbols('a:c', cls=Dummy)
```
But our labels usually mean something eg i,j will usually be occupied orbitals and a,b will be virtual orbitals. In the **particle-hole formalism (PHF)** we say occupied levels (particles) are below the Fermi surface and virtual ones (holes) above. We can impose these conditions on our labels by
```python
  i, j = symbols('i, j', below_fermi=True)
  a, b = symbols('a, b', above_fermi=True)
```
+ **creation and annihilation operators**\
To connect the occupied and virtual states we use creation operators (a<sup>&#8224;</sup>) and annihilation operators (a) (Shavitt & Bartlett 3.6).
```python
CreateFermion(p)*AnnihilateFermion(q)
or
Fd(p)*F(q)
```
+ An important concept for these operator is **normal order** (Shavitt & Bartlett 3.89), the normal order of these operators is one in which all the creation operators come to the left of the annihilation operators. We can achieve this by progressively swapping adjacent operators making a change of sign of the overall expression each time. Sympy has an operator to do this
```python
NO(AnnihilateFermion(p)*CreateFermion(q))
>>>-CreateFermion(q)*AnnihilateFermion(p)
```
+ Tensors. We define a **tensor**, say t<sup>ab</sup><sub>ij</sub> with sympy as
```python
v = AntiSymmetricTensor('v', (a, b), (i, j))
or
f = AntiSymmetricTensor('f', (a, ), (i, ))
```
+ **Normal Ordered Hamiltonian**\
We can now write code to form the normal ordered Hamiltonian (Shavitt & Bartlett 3.185)
```python
p, q, r, s = symbols('p:s', cls=Dummy)
i, j = symbols('i, j', below_fermi = True)
a, b = symbols('a, b', above_fermi= True)

#Fock tensor in normal order
f = AntiSymmetricTensor(operatorSymbol[0], (p, ), (q, ))
pq = NO(Fd(p)*F(q))

#2-electron repulsion tensor in normal order
g = AntiSymmetricTensor(operatorSymbol[4], (p, q), (r, s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))

#form Hamiltonian (normal order)
h = f * pq + Rational(1, 4) * g * pqsr
```
+ **Contractions**\
The contraction of operators A\*B is defined as A\*B - n\[A\*B] , where n\[&nbsp;] means the normal order (Shavitt & Bartlett 3.95) and is written ![image](https://user-images.githubusercontent.com/73105740/124299098-f3442f80-db54-11eb-9446-2b33005949b8.png). In sympy
```python
contraction(Fd(i),F(j))
```
+ **Wicks Theorem**\
The time-independent form of Wick’s theorem states: *A product of a string of creation and annihilation operators is equal to their normal product plus the sum of all possible normal products with contractions*. In sympy
```python
wicks(Fd(p)*F(q))
```
+ **Simplifying expressions**\
There are some useful functions which rationalise expressions. **evaluate_deltas(exp)**  does just that obeying Einstein summation convention. **substitute_dummies(exp)** this routine simplifys Add expressions containing terms which differ only due to dummy variables.

+ **Cluster Operators**\
The cluster operators are defined as T<sub>1</sub> = &Sigma; t<sup>a</sup><sub>i</sub> n\[a<sup>&#8224;</sup>i], T<sub>2</sub> = 0.25 &Sigma; t<sup>ab</sup><sub>ij</sub> n\[a<sup>&#8224;</sup>ib<sup>&#8224;</sup>j] etc (Shavitt & Bartlett 9.26-9.28, 9.29). In sympy we can do this as 
```python
    if 'S' in level:
        i = symbols('i', below_fermi=True, cls=Dummy)
        a = symbols('a', above_fermi=True, cls=Dummy)
        ts = AntiSymmetricTensor('ts', (a,), (i,))
        ai = NO(Fd(a)*F(i))
        t1 = ts * ai

    if 'D' in level:
        i, j = symbols('i,j', below_fermi=True, cls=Dummy)
        a, b = symbols('a,b', above_fermi=True, cls=Dummy)
        td = AntiSymmetricTensor('td', (a, b), (i, j))
        aibj = NO(Fd(a)*F(i)*Fd(b)*F(j))
        t2 = Rational(1, 4)*td*aibj

    if 'T' in level:
        i, j, k = symbols('i:k', below_fermi=True, cls=Dummy)
        a, b, c = symbols('a:c', above_fermi=True, cls=Dummy)
        tt = AntiSymmetricTensor('tt', (a, b, c), (i, j, k))
        aibjck = NO(Fd(a)*F(i)*Fd(b)*F(j)*Fd(c)*F(k))
        t3 = Rational(1, 36)*tt*aibjck

    if level == 'D':   return t2
    if level == 'SD':  return t1 + t2
    if level == 'SDT': return t1 + t2 + t3
``` 
+ **similarly-transformed Hamiltonian**\
This is defined as e<sup>-T</sup>**H**<sub>N</sub>e<sup>T</sup>, where **H**<sub>N</sub> is the normal-ordered Hamiltonian and T is a cluster operator. This is non-Hermitian giving rise to left and right eigenvectors.
+ **commutator**
The commutator bracket is defined \[A,B] = AB - BA and in sympy
```python
c = Commutator
```
+ **Permutations**
Sympy has a permutation operator 
```python
PermutationOperator(i,j)
```
+ **Baker-Campbell-Hausdorff**\
 The Baker-Campbell-Hausdorff expansion (BCH) gives an expansion of e<sup>-B</sup>Ae<sup>B</sup> (Shavitt & Bartlett 10.4) and applied to **H**<sub>N</sub> (Shavitt & Bartlett 10.5). A sympy implementation looks like
 ```python
 
    prettySymbols = {'above': 'defg','below': 'lmno', 'general': 'pqrst' }

    #commutator bracket
    c = Commutator
    
    bch = zeros(5)
    bch[0] = h
    for i in range(4):
        t  = clusterOperators(level)
        bch[i+1] = wicks(c(bch[i], t))
        bch[i+1] = evaluate_deltas(bch[i+1])          
        bch[i+1] = substitute_dummies(bch[i+1])

    #BCH series
    BCH = bch[0] + bch[1] + bch[2]/2 + bch[3]/6 + bch[4]/24 

    #tidy up and compact
    BCH = BCH.expand()
    BCH = evaluate_deltas(BCH)
    BCH = substitute_dummies(BCH, new_indices=True,
                                  pretty_indices=prettySymbols)
```
Here we have used a sympy array to hold the terms of the BCH expansion. Substitute_dummies just uses a user defined symbol dictionary to use in place of the internal dummies.
+ **Coupled-Cluster Doubles**
We'll see if we can use the above to get the CCD energy and amplitude expressions.
  + generate the normal-ordered Hamiltonian (see above)
  + do the BCH
  ```python
  cc = bakerCampbellHausdorff(h, 'D')
  ```
  + the cluster energy
  ```python
  w = wicks(cc, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
  symbol_rules = {'below':'ijklmno', 'above': 'abcdef', 'general':'pqrstu'}
  ccEnergy = substitute_dummies(w, new_indices=True, pretty_indices=symbol_rules)
  ```
  + the double amplitudes
  ```python
  expression = wicks(Fd(i)*Fd(j)*F(b)*F(a)*cc, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
  symbol_rules = {'below':'klmno', 'above': 'cdef', 'general':'pqrstu'}
  td = substitute_dummies(expression, new_indices=True, pretty_indices=symbol_rules)
  p = [PermutationOperator(i,j), PermutationOperator(a,b)]
  doubles = simplify_index_permutations(td, p)
  ```
+ The output\
  If we run the above code we get for
```
 **Energy**
 AntiSymmetricTensor(t, (_a, _b), (_i, _j))*AntiSymmetricTensor(v, (_i, _j), (_a, _b))/4
```

```
**Doubles Amplitudes** 
  AntiSymmetricTensor(f, (_k,), (i,))*AntiSymmetricTensor(t, (a, b), (j, _k))*PermutationOperator(i, j) 
- AntiSymmetricTensor(f, (a,), (_c,))*AntiSymmetricTensor(t, (b, _c), (i, j))*PermutationOperator(a, b) 
+ AntiSymmetricTensor(t, (_c, _d), (i, j))*AntiSymmetricTensor(t, (a, b), (_k, _l))*AntiSymmetricTensor(v, (_k, _l), (_c, _d))/4 
+ AntiSymmetricTensor(t, (_c, _d), (i, j))*AntiSymmetricTensor(v, (a, b), (_c, _d))/2 
+ AntiSymmetricTensor(t, (_c, _d), (j, _k))*AntiSymmetricTensor(t, (a, b), (i, _l))*AntiSymmetricTensor(v, (_k, _l), (_c, _d))*PermutationOperator(i, j)/2 
- AntiSymmetricTensor(t, (a, _c), (_k, _l))*AntiSymmetricTensor(t, (b, _d), (i, j))*AntiSymmetricTensor(v, (_k, _l), (_c, _d))*PermutationOperator(a, b)/2 
+ AntiSymmetricTensor(t, (a, _c), (i, _k))*AntiSymmetricTensor(t, (b, _d), (j, _l))*AntiSymmetricTensor(v, (_k, _l), (_c, _d))*PermutationOperator(i, j) 
+ AntiSymmetricTensor(t, (a, _c), (i, _k))*AntiSymmetricTensor(v, (b, _k), (j, _c))*PermutationOperator(a, b)*PermutationOperator(i, j) 
+ AntiSymmetricTensor(t, (a, b), (_k, _l))*AntiSymmetricTensor(v, (_k, _l), (i, j))/2 
+ AntiSymmetricTensor(v, (a, b), (i, j))
```
- - -
+ **Parsing the output**\
As you see the energy and amplitudes are each returned as a single string. The main tool we have is **args**, td.args returns a tuple of each term in the string. Each element in the tuple could be a term with multiple tensors (eg *+ AntiSymmetricTensor(t, (\_c,\_d), (i, j)) \* AntiSymmetricTensor(v, (a, b), (\_c, \_d))/2*) or a single tensor (eg *+ AntiSymmetricTensor(v, (a, b), (i, j))*). (It could in some situations be a constant). You can tell which sort it is using **isinstance**,eg isinstance(e, Mul), isinstance(e, AntiSymmetricTensor) or isinstance(e, Number).
```python
e = AntiSymmetricTensor(v, (a, b), (i, j))
isinstance(e, AntiSymmetricTensor)
>>>True
type(e)
>>>AntiSymmetricTensor
```
If we have a 'Mul' type, that is an expression with multiple elements, we can do another args on it. So say we do td.args\[0].args
```python
td.args[0]
>>>AntiSymmetricTensor(t, (a, b), (_k, _l))*AntiSymmetricTensor(v, (_k, _l), (i, j))/2
td.args[0].args
>>>(1/2, AntiSymmetricTensor(t, (a, b), (_k, _l)), AntiSymmetricTensor(v, (_k, _l), (i, j)))
```
This gives us another tuple (\[factor],Tensor,Tensor,...,Tensor,\[PermutationOperator]).\
If we have a AntiSymmetricTensor type we can do
```python
isinstance(td.args[0].args[1], AntiSymmetricTensor)
>>>True
td.args[0].args[1].symbol
>>>t
td.args[0].args[1].upper
>>>(a, b)
td.args[0].args[1].lower
>>>(_k, _l)
```
So we do an .args on the initial string, then do another args on all 'Mul' types. We now have either 'AntiSymmetricTensor', 'Number' or 'PermutationOperator' types. We have seen how to parse the AntiSymmetricTensor types, now 'PermutationOperator' types.
```
isinstance(td.args[4].args[4], PermutationOperator)
>>>True
td.args[4].args[4].args
>>>(i, j)
```
This is code to reduce the compound string to atoms (it might already be a 'Mul' in the case of energy)
```python
if isinstance(td, Mul):
    compound = [td]
else:
    compound = td.args
    
atoms = []
for element in compound:
    if isinstance(element, Mul):
        sub_atom = []
        for atom in element.args:
            sub_atom.append(atom)   
        atoms.append(sub_atom)
        
    elif isinstance(element, AntiSymmetricTensor):
        atoms.append([element])
        
    else:
        print('no current handler for this type: ',type(element))
        
```
- - -
+ Tests
  The output for CCD is, with Shavitt & Bartlett for comparison

 ![image](https://user-images.githubusercontent.com/73105740/124356517-51811900-dc0e-11eb-9bf9-6bdce09de87e.png)
![image](https://user-images.githubusercontent.com/73105740/124350850-ea07a100-dbee-11eb-82f1-b25d07f4d1ad.png)
 ![image](https://user-images.githubusercontent.com/73105740/124350750-6ea5ef80-dbee-11eb-9902-b8fa2ac7ff17.png)

These are the same allowing for permutation differences.

  The output for CCSD

 ![image](https://user-images.githubusercontent.com/73105740/124358259-6b732980-dc17-11eb-9064-c97d9dda0083.png)![image](https://user-images.githubusercontent.com/73105740/124358285-834aad80-dc17-11eb-89a0-97b71884cce7.png)
 
 Shavitt & Bartlett singles for comparison...
 
 ![image](https://user-images.githubusercontent.com/73105740/124358406-07049a00-dc18-11eb-84f9-36543c757ae9.png)
 
 Lets just check the last term in SB + v<sup>kl</sup><sub>cd</sub> t<sup>c</sup><sub>i</sub> t<sup>d</sup><sub>j</sub> t<sup>a</sup><sub>k</sub> t<sup>b</sup><sub>l</sub> matches our equation 22. First SB P(ij) v<sup>ab</sup><sub>cj</sub> t<sup>c</sup><sub>i</sub> is compared to our equation 8 which is -P(ij) v<sup>ab</sup><sub>jc</sub> t<sup>c</sup><sub>i</sub> or +P(ij) v<sup>ab</sup><sub>cj</sub> t<sup>c</sup><sub>i</sub> .

For CCSDT, the energy is the same. For singles we have an extra term 1/4 t<sup>**a**bc</sup><sub>**i**jk</sub> <jk||bc> compared to Shavitt & Bartlett 10.32 1/4 t<sup>**a**ef</sup><sub>**i**mn</sub> <mn||ef> where jkbc and mnef are repeated dummy indices.\
For the doubles there are 6 extra terms according to Shavitt & Bartlett 10.33, these are with our terms in {}, 
+ f<sub>me</sub> t<sup>abe</sup><sub>ijm</sub>&nbsp; &nbsp;{f<sub>kc</sub> t<sup>abc</sup><sub>ijk</sub> e<>c, m<>k}
+ 1/2P(ab) t<sup>aef</sup><sub>ijm</sub> <bm||ef>&nbsp; &nbsp;{1/2P(ab) t<sup>acd</sup><sub>ijk</sub> <bk||cd> e<>c. f<>d, m<>k}
+ -1/2P(ij) t<sup>abe</sup><sub>imn</sub> <mn||je>&nbsp; &nbsp;{-1/2P(ij) t<sup>abc</sup><sub>ikl</sub> <kl||jc> e<>c. m<>k, n<>l}
+ t<sup>e</sup><sub>m</sub> t<sup>fab</sup><sub>nij</sub> <mn||ef>&nbsp; &nbsp;{t<sup>c</sup><sub>k</sub> t<sup>abd</sup><sub>ijl</sub> <kl||cd> e<>c, m<>k, f<>d, n<>l}
+ -P(ab) t<sup>a</sup><sub>m</sub> t<sup>efb</sup><sub>inj</sub> <mn||ef>&nbsp; &nbsp; {1/2 P(ab) t<sup>a</sup><sub>k</sub> t<sup>bcd</sup><sub>ijl</sub> <kl||cd>  m<>k, n<>l, c<>e, d<>f}
+ -P(ij) t<sup>e</sup><sub>i</sub> t<sup>afb</sup><sub>mnj</sub> <mn||ef>&nbsp; &nbsp; {1/2 P(ij) t<sup>c</sup><sub>i</sub> t<sup>abd</sup><sub>jkl</sub> <kl||cd> m<>k, n<>l, c<>e, d<>f}

\*note symmetry of triples operator abc = -bac. The last 2 equations disagree by a factor 1/2 with Shavitt & Bartlett 10.33, but a check with [pdaggerq](https://github.com/edeprince3/pdaggerq/blob/master/examples/full_cc_codes/ccsdt.py) gives

+  1.0000 f(k,c)\*t3(c,a,b,i,j,k)
+	 0.5000 P(i,j)<l,k||c,j>\*t3(c,a,b,i,l,k)
+  0.5000 P(a,b)<k,a||c,d>\*t3(c,d,b,i,j,k)
+ -1.0000 <l,k||c,d>\*t1(c,k)\*t3(d,a,b,i,j,l)
+	-0.5000 P(i,j)<l,k||c,d>\*t1(c,j)\*t3(d,a,b,i,l,k) 
+	-0.5000 P(a,b)<l,k||c,d>\*t1(a,k)\*t3(c,d,b,i,j,l)

## code generation
As we've already parsed the expressions returned from sympy we could quite easily convert them into Python code. Here's a section from CCSDT triples...\
**note I define f<sup>a</sup><sub>i</sub> as f(i,a) - the subscript indices come first
and the superscript after so g<sup>ab</sup><sub>ij</sub> is f(i,j,a,b) this may be different from other sources.**
```python
    #  1.0000 * P(a,b) f(c,k) t1(k,a) t2(i,j,b,c) 
    t = 1.0000 * einsum('ck,ka,ijbc->ijab' ,f[v,o] ,t1 ,t2, optimize=True)
    T += t - t.swapaxes(3, 2)

    #  1.0000 * P(i,j) f(c,k) t1(i,c) t2(j,k,a,b) 
    t = 1.0000 * einsum('ck,ic,jkab->ijab' ,f[v,o] ,t1 ,t2, optimize=True)
    T += t - t.swapaxes(1, 0)

    #  1.0000 * P(a,b)P(i,j) <i,c||a,k> t2(j,k,b,c) 
    t = 1.0000 * einsum('icak,jkbc->ijab' ,g[o,v,v,o] ,t2, optimize=True)
    T += t - t.swapaxes(3, 2) - t.swapaxes(1, 0) + t.swapaxes(3, 2).swapaxes(1, 0)
```
pdaggerq uses einsum to do the permutations, its neat and shows whats happening clearly but I wonder if it's as efficient as swapaxes (or indeed transpose).
The code generator also writes a subroutine header and return statement.
```python
def cogusDouble(f, g, o, v, t1=None, t2=None, t3=None):

    '''
        COGUS generated level [SDT] on 07/08/21   
    '''
    from numpy import einsum, swapaxes
```
which is the same for all routines, just the fock and 2-electron repulsion integrals followed by occupied and virtual slices, then the amplitudes if they exist.

## test routine
I've written a routine to test the generated cluster amplitude. This
+ takes a **type**- C (for cluster), **level**- D/SD/SDT, molecule- h2o/ch4/acetaldehyde and basis- sto-3g/6-31g/dz as command line arguments. Valid combinations are

|  level | molecule  | basis |
|--------|-----------|-------|
|   D    |   h2o     | sto-3g|
| SD     | h2o       | sto-3g/6-31g/dz |
| SD     | acetaldehyde/ch4 | sto-3g |
| SDT    |  h2o   | sto-3g  |

Tested against pdaggerq code where available ie SD and SDT. 
+ checks molecule and basis against those defined in the harpy project file (first entry)
+ reads mints file for molecule and basis, reads .cdf file of python code for cluster level (eg cogus_C_SD.cdf)
+ concatenates cluster data read in, with internal code to run loop over amplitude generation.
+ exec's code
+ compares energy correction with reference values.

So no cluster code is contained within the test program as it gets it from an external file previously generated by cogus.
If cogus was fast (as I suspect pdaggerq is) then the cluster code could be generated on-the-fly rather than being read in from file. Anyway the cogus generated code for D, SD and SDT amplitudes agree's with reference values to better than 1e-8.

## Linear Coupled Cluster
It is easy to implement LCCD and LCCSD by simply truncting the Baker-Campbell-Hausdorff expansion after 2 terms (ie only allowing linear contributions). To allow for this I've defined a new type 'L' to go along with 'C'. Python linear coupled-cluster code files should be named 'cogus_L_SD.cdf' for example. cogusValidate has been updated to check LCCD abd LCCSD on H<sub>2</sub>O in STO-3G basis. To implement this 'type' has been added to the arguments of bakerCampbellHausdorff then a dictionary {'C':4,'L':1} is used to determine the number of loops (terms) of the expansion to include.

## CC2
This has now been implemented as type 'A' with level '2', called as *cogusGenerate.py A 2 E+S+D True*. CC2 has identical E and S code to CCSD however the doubles is defined as <0| i<sup>+</sup>j<sup>+</sup>ab (e<sup>-t1-t2</sup>**f**e<sup>t1+t2</sup> + e<sup>-t1</sup>**v**e<sup>t1</sup>) |0> . To implement this an extra option is returned by the cluster routine 'if level == 'S':   return T1'. The doubles evaluation now become
```python
            #Cluster doubles amplitude
            if type in ['A']: cc = bakerCampbellHausdorff(f, type, 'D') + bakerCampbellHausdorff(v, type, 'S')
            w = wicks(Fd(i)*Fd(j)*F(b)*F(a)*cc, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
            symbolRules = {'below':'klmno', 'above': 'cdef', 'general':'pqrstu'}
            td = substitute_dummies(w, new_indices=True, pretty_indices=symbolRules)
            p = [PermutationOperator(i,j), PermutationOperator(a,b)]

            doubles = simplify_index_permutations(td, p)
```
cogusValidate now has options 'A 2 h2o sto-3g' and 'A 2 h2o dz' with validation against Psi4.

## CC3
This has now been implemented as type 'A' with level '3', called as *cogusGenerate.py A 3 E+S+D+T True*. CC3 has identical E, S and D code to CCSDT however the triples is defined as <0| i<sup>+</sup>j<sup>+</sup>k<sup>+</sup>abc (e<sup>-t1-t2-t3</sup>**f**e<sup>t1+t2+t3</sup> + i<sup>+</sup>a e<sup>-t1</sup>**v**e<sup>t1</sup> + \[v, t2] + \[\[v, t1], t2] + (1/2)\[\[\[v, t1], t1], t2] + (1/6)\[\[\[\[v, t1], t1], t1], t2]) |0> . This is implemented in triples as
```python
if (type == 'A') and (level == 'SDT'):
   cc = bakerCampbellHausdorff(f, 'C', 'SDT') + bakerCampbellHausdorff(v, 'C', 'S') + v +  \
        Commutator(v, clusterOperators('D')) + Commutator(Commutator(v,clusterOperators('S')),\
        clusterOperators('D')) + Rational(1,2)*Commutator(Commutator(Commutator(v, clusterOperators('S')),\
	clusterOperators('S')),clusterOperators('D')) + Rational(1,6)*Commutator(Commutator(Commutator(Commutator(v, \
	clusterOperators('S')), clusterOperators('S')),clusterOperators('S')),clusterOperators('D'))
```
cogusValidate now has options 'A 3 h2o sto-3g' and 'A 2 h2o dz' with validation against Psi4.

## CCSD(T) 
This has been implemented as type 'C' with level 'SDt'. The S and D amplitudes are the same as CCSD but the triples amplitudes
are given by <0| i<sup>+</sup>j<sup>+</sup>k<sup>+</sup>abc (e<sup>-t3</sup>**f**e<sup>t3</sup> + \[v, t2]) |0> this is added non-iteratively. The perturbation energy correction is then 
(l1+l2) operating on \[v, t3] where l1 and l2 are the Lagrange multiplier operators (transpose of t amplitudes). This is implemented for triples and energy as... 
```python
         cc = bakerCampbellHausdorff(f, 'C', 'T') + Commutator(v, clusterOperators('D'))
         cc = wicks(Fd(i)*Fd(j)*Fd(k)*F(c)*F(b)*F(a)*cc ,keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
         symbolRules = {'below':'lmno', 'above': 'def', 'general':'pqrstu'}

         cc = substitute_dummies(cce, new_indices=True, pretty_indices=symbolRules)
         p = [PermutationOperator(i,j), PermutationOperator(a,b),PermutationOperator(j,k), PermutationOperator(b,c)]  
         mixture['T'] = simplify_index_permutations(cc, p)

        #singles 
         i = symbols('i', below_fermi=True, cls=Dummy)
         a = symbols('a', above_fermi=True, cls=Dummy)
         ls = AntiSymmetricTensor('l1', (i,), (a,))
         ai = NO(Fd(i)*F(a))
         l1 = ls * ai
	 
         #doubles
         i, j = symbols('i,j', below_fermi=True, cls=Dummy)
         a, b = symbols('a,b', above_fermi=True, cls=Dummy)
         ld = AntiSymmetricTensor('l2', (i, j), (a, b))
         aibj = NO(Fd(i)*F(a)*Fd(j)*F(b))
         l2 = Rational(1, 4)*ld*aibj

         cc =   Commutator(v, clusterOperators('T'))
         lop = l1 + l2
         w = wicks(lop*cc, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
         symbolRules = {'below':'ijklmno', 'above': 'abcdef', 'general':'pqrstu'}

         mixture['E'] = substitute_dummies(w, new_indices=True, pretty_indices=symbolRules)

```
cogusValidate now has options 'C SDt h2o sto-3g', 'h2o 6-31g', 'ch4 sto-3g'

## CCSD Response Density Matrices
Shavitt & Bartlett 11.88 defines the response density matrix as (&gamma;)<sub>qp</sub> = <0| (1 + &Lambda;)e<sup>-T</sup>{p<sup>+</sup>q}e<sup>T</sup> |0>\
here &Lambda; are de-excitation operators defined in a similar way to T ie &Lambda;<sub>n</sub>=1/(n!)<sup>2</sup> &Sigma;&lambda;<sup>i..n</sup><sub>a..f</sub> (i<sup>+</sup>a..n<sup>+</sup>f}.\
We can generalise the density matrix equation to include EOM cases
|   method  |   &Lambda;<sub>0</sub>   |  &Lambda;<sub>1</sub> |  &Lambda;<sub>2</sub> |
|-----------|--------------------------|-----------------------|-----------------------|
|    CC     |                1         | &lambda;<sup>i</sup><sub>a</sub> i<sup>+</sup>a | (1/4)&lambda;<sup>jk</sup><sub>bc</sub> j<sup>+</sup>k<sup>+</sup>cb |
|    EE     |                0         | &lambda;<sup>i</sup><sub>a</sub> i<sup>+</sup>a | (1/4)&lambda;<sup>jk</sup><sub>bc</sub> j<sup>+</sup>k<sup>+</sup>cb |
|    IP     |                0         | &lambda;<sup>i</sup> i<sup>+</sup> | (1/2)&lambda;<sup>jk</sup><sub>a</sub> j<sup>+</sup>k<sup>+</sup>a |
|    EA     |                0         | &lambda;<sub>a</sub> a | (1/2)&lambda;<sup>i</sup><sub>bc</sub> i<sup>+</sup>cb |

The basic EE form has i<sup>+</sup>a. if we go to IP we lose the hole ie i<sup>+</sup> or if we go to electron affinity we lose the particle a. There are also double IP (DIP) and double EA (DEA).\
In EOM the density matrices are given by <0| L<sub>k</sub> e<sup>-T</sup> {p<sup>+</sup>q} e<sup>T</sup> R<sub>l</sub> |0> Shavitt and Bartlett 13.30. Here L and R are the left and right eigenvectors Shavitt & Bartlett 13.14 and 13.9, they represent de-excitations and excitations, respectively. They form of bi-orthonormal set Shavitt & Bartlett 13.16 - more of EOM later.

We can implement as 
```python
method = 'CC'

def de_excitations(method):

   i, j, k = symbols('i:k' ,below_fermi=True, cls=Dummy)
   a ,b, c = symbols('a:c' ,above_fermi=True, cls=Dummy)   

   if method == 'IP':
       return [0, AntiSymmetricTensor('l',(i,),())*Fd(i), Rational(1, 2)* \
                  AntiSymmetricTensor('l',(j,k),(a,))*Fd(j)*Fd(k)*F(a)]  
   elif method == 'EA':
       return [0, AntiSymmetricTensor('l',(),(a,))*F(a), Rational(1, 2)* \
                  AntiSymmetricTensor('l',(i,),(b,c))*Fd(i)*F(c)*F(b)]
   elif method == 'EE':
       return [0, AntiSymmetricTensor('l',(i,),(a,))*Fd(i)*F(a), Rational(1, 4)* \
                  AntiSymmetricTensor('l',(j,k),(b,c))*Fd(j)*Fd(k)*F(c)*F(b)]
   elif method == 'CC':
       return [1, AntiSymmetricTensor('l',(i,),(a,))*Fd(i)*F(a), Rational(1, 4)* \
                  AntiSymmetricTensor('l',(j,k),(b,c))*Fd(j)*Fd(k)*F(c)*F(b)]
  
L = sum(de_excitations(method))

#as example D(oo)
cc = bakerCampbellHausdoff(Fd(i)*F(j),'C','SD')
evaluate_deltas(wicks(L*cc, keep_only_fully_contracted = True))
cc = substitute_dummies(cc, new_indices=True, pretty_indices = {'below':  'klmno','above':  'abcde'})

#and D(vv)
cc = bakerCampbellHausdoff(Fd(a)*F(b),'C','SD')
evaluate_deltas(wicks(L*cc, keep_only_fully_contracted = True))
cc = substitute_dummies(cc, new_indices=True, pretty_indices = {'below':  'ijklmn','above':  'cdefg'})
```
This results in
```
oo=  KroneckerDelta(i, j) - AntiSymmetricTensor(l, (j,), (_a,))*AntiSymmetricTensor(t, (_a,), (i,)) \
   - AntiSymmetricTensor(l, (j, _k), (_a, _b))*AntiSymmetricTensor(t, (_a, _b), (i, _k))/2

vv=  AntiSymmetricTensor(l, (_i,), (a,))*AntiSymmetricTensor(t, (b,), (_i,)) + \
     AntiSymmetricTensor(l, (_i, _j), (a, _c))*AntiSymmetricTensor(t, (b, _c), (_i, _j))/2
```
or &delta;<sub>ij</sub> - l<sup>j</sup><sub>a</sub> t<sup>a</sup><sub>i</sub> - 0.5 l<sup>jk</sup><sub>ab</sub> t<sup>ab</sup><sub>ik</sub>\
and l<sup>i</sup><sub>a</sub> t<sup>b</sup><sub>i</sub> + 0.5 l<sup>ij</sup><sub>ac</sub> t<sup>bc</sup><sub>ij</sup>
