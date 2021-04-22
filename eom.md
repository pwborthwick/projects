# EOM_CCSD
The paper [Simplified methods for equation-of-motion coupled-cluster excited state calculations - Steven R. Gwaltney, Marcel Nooijen, Rodney J. Bartlett](https://notendur.hi.is/agust/rannsoknir/papers/cpl248-189-96.pdf) gives the following equations for the partitioning of the EOM-CCSD hamiltonian 

H<sub>SS</sub> = &Sigma; F<sub>ac</sub> - &Sigma; F<sub>ik</sub> + &Sigma; W<sub>akic</sub>

H<sub>SD</sub> = &Sigma; F<sub>ld</sub> + 0.5 &Sigma; W<sub>alcd</sub> - 0.5 &Sigma; W<sub>klid</sub>

H<sub>DS</sub> =*P*(ab) W<sub>kaij</sub> + *P*(ij) W<sub>abej</sub> + *P*(ab) W<sub>bmfe</sub> t<sup>af</sup><sub>ij</sub> - *P*(ij) W<sub>mnje</sub> t<sup>ab</sup><sub>in</sub> 

- - - -
#### H<sub>SS</sub>
Equations for terms taken from  [J. Chem. Phys. 98, 7029 (1993); https://doi.org/10.1063/1.46474698, 7029© 1993 American Institute of Physics.The equation of motion coupled-clustermethod. A systematic biorthogonal approach to molecular excitation energies, transition probabilities, and excited state properties](https://www.theochem.ru.nl/files/local/jcp-98-7029-1993.pdf) 

+ F<sub>**ac**</sub> = *f*<sub>**ac**</sub> &delta;<sub>**ik**</sub> - t<sub>**a**</sub><sup>m</sup> *f*<sub>m**c**</sub>&delta;<sub>**ik**</sub> + t<sup>em</sup>g<sub>m**a**e**c**</sub>&delta;<sub>**ik**</sub> - 0.5 t<sup>e</sup><sub>**a**</sub><sup>nm</sup>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub> - t<sup>e</sup><sub>n</sub>t<sup>**a**</sup><sub>m</sub>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub>

    + [1]  +*f*<sub>**ac**</sub> &delta;<sub>**ik**</sub>
    + [2]  -t<sub>**a**</sub><sup>m</sup> *f*<sub>m**c**</sub>&delta;<sub>**ik**</sub>
    + [3]  +t<sup>em</sup>g<sub>m**a**e**c**</sub>&delta;<sub>**ik**</sub> 
    + [4]  -0.5 t<sup>e</sup><sub>**a**</sub><sup>nm</sup>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub> 
    + [5]  -t<sup>e</sup><sub>n</sub>t<sup>**a**</sup><sub>m</sub>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub>
   
+ F<sub>**ik**</sub> = *f*<sub>**ik**</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>**k**</sub> *f*<sub>**i**e</sub>&delta;<sub>**ac**</sub> + t<sup>em</sup>g<sub>**i**m**k**e</sub>&delta;<sub>**ac**</sub> + 0.5 t<sup>ef</sup><sub>**k**</sub><sup>m</sup>g<sub>**i**mef</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>m</sub>t<sup>f</sup><sub>**k**</sub>g<sub>**i**mef</sub>&delta;<sub>**ac**</sub> 

    + [6]  +*f*<sub>**ik**</sub>&delta;<sub>**ac**</sub>
    + [7]  +t<sup>e</sup><sub>**k**</sub> *f*<sub>**i**e</sub>&delta;<sub>**ac**</sub>
    + [8]  +t<sup>em</sup>g<sub>**i**m**k**e</sub>&delta;<sub>**ac**</sub> 
    + [9]  +0.5 t<sup>ef</sup><sub>**k**</sub><sup>m</sup>g<sub>**i**mef</sub>&delta;<sub>**ac**</sub>
    + [10] +t<sup>e</sup><sub>m</sub>t<sup>f</sup><sub>**k**</sub>g<sub>**i**mef</sub>&delta;<sub>**ac**</sub> 

+ W<sub>**akic**</sub> = g<sub>**akic**</sub> + t<sup>e</sup><sub>**i**</sub>g<sub>**ak**e**c**</sub> - t<sup>**a**</sup><sub>m</sub>g<sub>m**ki**c</sub> - (t<sup>e**a**</sup><sub>**i**m</sub> + t<sup>e</sup><sub>**i**</sub>t<sup>**a**</sup><sub>m</sub>)g<sub>m**k**e**c**</sub>

    + [11] +g<sub>**akic**</sub>
    + [12] +t<sup>e</sup><sub>**i**</sub>g<sub>**ak**e**c**</sub>
    + [13] -t<sup>**a**</sup><sub>m</sub>g<sub>m**ki**c</sub>
    + [14] -t<sup>e**a**</sup><sub>**i**m</sub>g<sub>m**k**e**c**</sub>
    + [15] -t<sup>e</sup><sub>**i**</sub>t<sup>**a**</sup><sub>m</sub>g<sub>m**k**e**c**</sub>
 
- - -
#### H<sub>SD</sub>
+ F<sub>ld</sub> = f<sub>**ld**</sub>&delta;<sub>**ik**</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>m</sub>g<sub>**k**m**c**e</sub>&delta;<sub>**il**</sub>&delta;<sub>**ad**</sub>

    + [16] +f<sub>**ld**</sub>&delta;<sub>**ik**</sub>&delta;<sub>**ac**</sub>
    + [17] +t<sup>e</sup><sub>m</sub>g<sub>**k**m**c**e</sub>&delta;<sub>**il**</sub>&delta;<sub>**ad**</sub>
    
+  W<sub>**alcd**</sub> = g<sub>**alcd**</sub>&delta;<sub>**ik**</sub> - t<sup>**a**</sup><sub>m</sub>g<sub>m**lcd**</sub>&delta;<sub>**ik**</sub>

    + [18] +g<sub>**alcd**</sub>&delta;<sub>**ik**</sub> 
    + [19] -t<sup>**a**</sup><sub>m</sub>g<sub>m**lcd**</sub>&delta;<sub>**ik**</sub>

+ W<sub>klid</sub> = g<sub>**klid**</sub> </sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>**i**</sub>g<sub>**kl**e**d**</sub></sub>&delta;<sub>**ac**</sub>

    + [20] +g<sub>**klid**</sub> </sub>&delta;<sub>**ac**</sub>
    + [21] +t<sup>e</sup><sub>**i**</sub>g<sub>**kl**e**d**</sub></sub>&delta;<sub>**ac**</sub>
    
*There is disagreement between reference [2] and [Coupled-cluster calculations of nuclear magnetic resonance chemical shifts](www2.chemia.uj.edu.pl/~migda/Literatura/pdf/JCP03561.pdf) we have taken reference [3] which agrees with coding in psi4numpy/pyscf. Reference [2] has g<sub>likd</sub> + t<sup>d</sup><sub>k</sub>g<sub>lied</sub> and reference [3] g<sub>klid</sub> + t<sup>e</sup><sub>i</sub>g<sub>kled</sub>*

- - -
### H<sub>DS</sub>
+ W<sub>kaij</sub> = g<sub>**kaij**</sub>&delta;<sub>**bc**</sub> + *P*(ij) t<sup>e**a**</sup><sub>m**j**</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**bc**</sub> + 0.5&tau;<sup>ef</sup><sub>**ij**</sub>g<sub>**ka**ef</sub>&delta;<sub>**bc**</sub> - t<sup>**a**</sup><sub>m</sub>W<sub>**k**m**ij**</sub>&delta;<sub>**bc**</sub> + +*P*(ij) t<sup>e</sup><sub>**j**</sub> (g<sub>**ia**e**k**</sub> - t<sup>**a**f</sup><sub>m**k**</sub>g<sub>**i**mef</sub>)&delta;<sub>**bc**</sub> +t<sup>ea</sup><sub>ij</sub>F<sub>ke</sub>&delta;<sub>**bc**</sub>

    + [22] +g<sub>**kaij**</sub>&delta;<sub>**bc**</sub>
    + [23] +*P*(ij) t<sup>e**a**</sup><sub>m**j**</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**bc**</sub>
    + [--] +0.5&tau;<sup>ef</sup><sub>**ij**</sub>g<sub>**ka**ef</sub>&delta;<sub>**bc**</sub> 
        + [24] +0.5t<sup>ef</sup><sub>**ij**</sub>g<sub>**ka**ef</sub>&delta;<sub>**bc**</sub> 
        + [25] +t<sup>e</sup><sub>**i**</sub>t<sup>f</sup><sub>**j**</sub>g<sub>**ka**ef</sub>&delta;<sub>**ac**</sub>
    + [--] -t<sup>a</sup><sub>m</sub>W<sub>kmij</sub>&delta;<sub>**bc**</sub>
        + W<sub>kmij</sub> =  g<sub>kmij</sub> + *P*(ij) t<sup>e</sup><sub>j</sub>g<sub>kmie</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub>g<sub>kmef</sub> 
            + [26] -t<sup>a</sup><sub>m</sub>g<sub>kmij</sub> &delta;<**bc**</sub>
            + [27] +t<sup>**a**</sup><sub>m</sub> *P*(ij) t<sup>e</sup><sub>**j**</sub>g<sub>**k**me**i**</sub> &delta;<**bc**</sub>
            + [--] -0.5&tau;<sup>**a**</sup><sub>m</sub>t<sup>ef</sup><sub>ij</sub>g<sub>kmef</sub> &delta;<sub>**bc**</sub>
                + [28] -0.5t<sup>**a**</sup><sub>m</sub>t<sup>ef</sup><sub>ij</sub>g<sub>kmef</sub> &delta;<sub>**bc**</sub>
                + [29] -t<sup>**a**</sup><sub>m</sub>t<sup>e</sup><sub>**i**</sub>t<sup>f</sup><sub>**j**</sub>g<sub>**k**mef</sub> &delta;<sub>**bc**</sub>
    + [--] +*P*(ij) t<sup>e</sup><sub>**j**</sub> (g<sub>**ia**e**k**</sub> - t<sup>**a**f</sup><sub>m**k**</sub>g<sub>**i**mef</sub>)&delta;<sub>**bc**</sub>
        + [30] +*P*(ij) t<sup>e</sup><sub>**i**</sub>g<sub>**ka**e**j**</sub>&delta;<sub>**bc**</sub> 
        + [31] -*P*(ij) t<sup>e</sup><sub>**i**</sub> t<sup>**a**f</sup><sub>m**j**</sub>g<sub>**k**mef</sub>&delta;<sub>**bc**</sub> 
    + [32] +t<sup>ea</sup><sub>ij</sub>F<sub>ke</sub>&delta;<sub>**bc**</sub>



### H<sub>DS</sub>

From reference [2] table I 
W<sub>iajk</sub> = g<sub>iajk</sub> + *P*(ij)t<sup>ae</sup><sub>km</sub>g<sub>imje</sub> + 0.5&tau;<sup>ef</sup><sub>jk</sub>g<sub>iaef</sub> + t<sup>a</sup><sub>m</sub>W<sub>imjk</sub>  - *P*(ij) t<sup>e</sup><sub>j</sub> (g<sub>iaek</sub> - t<sup>af</sup><sub>mk</sub>g<sub>imef</sub> ) + t<sup>ae</sup><sub>jk</sub>F<sub>ie</sub>

            
      
Again the sign difference between references \[2] and \[3].
psi4numpy has the code
```python
def build_Wmbij():
    """Eqn 9 from [Gauss:1995:3561], Table III (b)"""
    Wmbij = ccsd.get_MO('ovoo').copy()
    
    Wmbij -= np.einsum('me,ijbe->mbij', Fme, t2)
    Wmbij -= np.einsum('nb,mnij->mbij', t1, Wmnij)
   
    temp_tau = ccsd.build_tau() 
    Wmbij += 0.5 * np.einsum('mbef,ijef->mbij', ccsd.get_MO('ovvv'),temp_tau)
   
    Pij = np.einsum('jnbe,mnie->mbij', t2, ccsd.get_MO('ooov'))
    Wmbij += Pij
    Wmbij -= Pij.swapaxes(2, 3)

    temp_mbej = ccsd.get_MO('ovvo').copy()
    temp_mbej -= np.einsum('njbf,mnef->mbej', t2, ccsd.get_MO('oovv'))
    Pij = np.einsum('ie,mbej->mbij', t1, temp_mbej)
    Wmbij += Pij
    Wmbij -= Pij.swapaxes(2, 3)
    return Wmbij
```
From this we have W<sub>mbij</sub> = g<sub>mbij</sub> **-** t<sup>be</sup><sub>ij</sub>F<sub>me</sub> **-** t<sup>b</sup><sub>n</sub>W<sub>mnij</sub> **+** 0.5&tau;<sup>ef</sup><sub>ij</sub>g<sub>mbef</sub> **+** *P*(ij) t<sup>be</sup><sub>jn</sub>g<sub>mnie</sub> **+** *P*(ij) t<sup>e</sup><sub>i</sub>(g<sub>mbej</sub> **-** t<sup>bf</sup><sub>nj</sub>g<sub>mnef</sub>

Hence the second *P*(ij) term is with a plus sign

The sign on the W<sub>imjk</sub> term is given by psi4numpy and reference [3] as **-**.

                                     

 From [3] (and psi4numpy) W<sub>mbij</sub> we have -t<sup>be</sup><sub>ij</sub>F<sub>me</sub> so for W<sub>maij</sub> we have -t<sup>ae</sup><sub>ij</sub>F<sub>me</sub>
re-index as -t<sup>**a**e</sup><sub>**ij**</sub>F<sub>**k**e</sub> and swapping index on t to +t<sup>ea</sup><sub>ij</sub>F<sub>ke</sub>&delta;<sub>**bc**</sub>
                  
                      HDS[iajb,kc] += (b==c)*fs[k,e]*td[e,a,i,j] - \
                                      (a==c)*fs[k,e]*td[e,b,i,j]  # 29
                                      
 Now for *P*(ij) W<sub>abej</sub>  terms, reference [3] gives W<sub>abei</sub> = g<sub>abei</sub> - *P*(ab) g<sub>mbef</sub>t<sup>af</sup><sub>mi</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub>g<sub>mnei</sub> + t<sup>f</sup><sub>i</sub>W<sub>abef</sub> - *P*(ab) t<sup>a</sup><sub>m</sub>(g<sub>mbei</sub> - t<sup>bf</sup><sub>ni</sub>g<sub>mnef</sub>
 
g<sub>abei</sub> re-index (e->c)   gives g<sub>**abci**</sub>&delta;<sub>**ik**</sub>       

                      HDS[iajb,kc] += (i==k)*ints[a,b,c,j] - \
                                      (j==k)*ints[a,b,c,i] # 22
                                      
-*P*(ab) g<sub>mbef</sub>t<sup>af</sup><sub>mi</sub> re-indexing (e->c,f->e) -g<sub>mbce</sub>t<sup>ae</sup><sub>mi</sub> interchanging indices gives

+g<sub>**b**m**c**e</sub>t<sup>**a**e</sup><sub>m**i**</sub>&delta;<sub>**jk**</sub>

                       (j==k)*ints[b,m,c,e]*td[e,a,m,i] # 30
                     - (j==k)*ints[a,m,c,e]*td[e,b,m,i] 

                       HDS[iajb,kc] += -(i==k)*ints[b,m,c,e]*td[e,a,m,j] + \
                                        (i==k)*ints[a,m,c,e]*td[e,b,m,j] - \
                                    
+0.5&tau;<sup>ab</sup><sub>mn</sub>g<sub>mnei</sub> = 0.5t<sup>ab</sup><sub>mn</sub>g<sub>mnei</sub> + t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>g<sub>mnei</sub>

0.5t<sup>ab</sup><sub>mn</sub>g<sub>mnei</sub> re-index 0.5t<sup>ab</sup><sub>mn</sub>g<sub>mncj</sub>
   
                      HDS[iajb,kc] += 0.5*(i==k)*ints[m,n,c,j]*td[a,b,m,n] - \
                                      0.5*(j==k)*ints[m,n,c,i]*td[a,b,m,n]  # 33 

t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>g<sub>mnei</sub> re-index t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>g<sub>mncj</sub>
  
                      HDS[iajb,kc] += (i==k)*ts[a,m]*ts[b,n]*ints[m,n,c,j] - \
                                      (j==k)*ts[a,m]*ts[b,n]*ints[m,n,c,i]  # 37

+t<sup>f</sup><sub>j</sub>W<sub>abef</sub>    from reference \[2] W<sub>abcd</sub> = g<sub>abcd</sub> - *P*(ab) t<sup>b</sup><sub>m</sub>g<sub>amcd</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub>g<sub>mncd</sub>

+t<sup>f</sup><sub>j</sub>g<sub>abef</sub> is then +t<sup>f</sup><sub>j</sub>g<sub>abef</sub> then +t<sup>e</sup><sub>j</sub>g<sub>abce</sub>

                      HDS[iajb,kc] += (i==k)*ints[a,b,c,e]*ts[e,j] - \
                                      (j==k)*ints[a,b,c,e]*ts[e,i] # 24

                                      

