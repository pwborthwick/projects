## Cluster Operator Generation Using Sympy

Sympy has a module called secondquant and I want to use this to automatically generate some cluster operators.
+ Generate some symbols\
  Quantum chemistry is pretty heavy on labels, we have lots of tensors, so how do using sympy we define some labels to use with our tensors. Unlike in Python where we just just use the variables we want without declaration, in sympy we must declare our variables/labels
  
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
+ An important concept of these operator is **normal order** (Shavitt & Bartlett 3.89), the normal order of these operators is one in which all the creation operators come to the left of the annihilation operators. We can achieve this by progressively swapping adjacent operators making a change of sign of the overall expression each time. Sympy has an operator to do this
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
+ Simplifying expreesions\
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
Simpy has a permutation operator 
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
   = simplify_index_permutations(td, p)
  ```
+ The output\
  If we run the above code we get for\
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
As you see the energy and amplitudes are returned as a single string. The main tool we have is **args**, td.args returns a tuple of each term in the string. Each element in the tuple could be a term with multiple tensors (eg *+ AntiSymmetricTensor(t, (\_c,\_d), (i, j)) \* AntiSymmetricTensor(v, (a, b), (\_c, \_d))/2*) or a single tensor (eg *+ AntiSymmetricTensor(v, (a, b), (i, j))*). (It could in some situations be a constant). You can tell which sort it is using **isinstance**,eg isinstance(e, Mul), isinstance(e, AntiSymmetricTensor) or isinstance(e, Number).
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




            
