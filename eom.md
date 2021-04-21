**# EOM-CCSD**
The paper [Simplified methods for equation-of-motion coupled-cluster excited state calculations - Steven R. Gwaltney, Marcel Nooijen, Rodney J. Bartlett](https://notendur.hi.is/agust/rannsoknir/papers/cpl248-189-96.pdf) gives the following equations for the partitioning of the EOM-CCSD hamiltonian 

H<sub>SS</sub> = &Sigma; F<sub>ae</sub> - &Sigma; F<sub>mi</sub> + &Sigma; W<sub>amie</sub>

H<sub>SD</sub> = &Sigma; F<sub>me</sub> + 0.5 &Sigma; W<sub>amef</sub> - 0.5 &Sigma; W<sub>mnie</sub>

H<sub>DS</sub> =*P*(ab) W<sub>maij</sub> + *P*(ij) W<sub>abej</sub> + *P*(ab) W<sub>bmfe</sub> t<sup>af</sup><sub>ij</sub> - *P*(ij) W<sub>mnje</sub> t<sup>ab</sup><sub>in</sub> 

Equations for terms taken from  [J. Chem. Phys. 98, 7029 (1993); https://doi.org/10.1063/1.46474698, 7029© 1993 American Institute of Physics.The equation of motion coupled-clustermethod. A systematic biorthogonal approach to molecular excitation energies, transition probabilities, and excited state properties](https://www.theochem.ru.nl/files/local/jcp-98-7029-1993.pdf) 

### H<sub>SS</sub>
From reference [2] table I

F<sub>ab</sub> = *f*<sub>ab</sub> - t<sub>a</sub><sup>m</sup> *f*<sub>mb</sub> + t<sup>fm</sup>g<sub>mafb</sub> - 0.5 &tau; <sup>e</sup><sub>a</sub><sup>nm</sup>g<sub>mnbe</sup>

We can re-index (b->c, f->e) as F<sub>**ac**</sub> = *f*<sub>**ac**</sub> &delta;<sub>**ik**</sub> - t<sub>**a**</sub><sup>m</sup> *f*<sub>m**c**</sub>&delta;<sub>**ik**</sub> + t<sup>em</sup>g<sub>m**a**e**c**</sub>&delta;<sub>**ik**</sub> - 0.5 t<sup>e</sup><sub>**a**</sub><sup>nm</sup>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub> - t<sup>e</sup><sub>n</sub>t<sup>**a**</sup><sub>m</sub>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub>

    HSS[ia,kc] += fs[<a>,<c>]*(i==k)                        # 1 
    HSS[ia,kc] += -(i==k)*fs[m,<c>]*ts[<a>,m]               # 9
    HSS[ia,kc] += ints[<a>,m,<c>,e]*ts[e,m]*(i==k)          # 4 
    HSS[ia,kc] += -0.5*(i==k)*ints[m,n,<c>,e]*td[<a>,e,m,n] # 10
    HSS[ia,kc] += -(i==k)*ints[m,n,<c>,f]*ts[<a>,m]*ts[f,n] # 15

From reference [2] table I

F<sub>ji</sub> = *f*<sub>ji</sub> + t<sup>e</sup><sub>i</sub> *f*<sub>je</sub> + t<sup>em</sup>g<sub>jmie</sub> + 0.5 &tau; <sup>efm</sup><sub>i</sub>g<sub>jmef</sub>

We can re-index (j->i, i->k) as F<sub>**ik**</sub> = *f*<sub>**ik**</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>**k**</sub> *f*<sub>**i**e</sub>&delta;<sub>**ac**</sub> + t<sup>em</sup>g<sub>**i**m**k**e</sub>&delta;<sub>**ac**</sub> + 0.5 t<sup>ef</sup><sub>**k**</sub><sup>m</sup>g<sub>**i**mef</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>m</sub>t<sup>f</sup><sub>**k**</sub>g<sub>**i**mef</sub>&delta;<sub>**ac**</sub> 



    HSS[ia,kc] += -fs[<i>,<k>]*(a==c)                       # 2
    HSS[ia,kc] += -(a==c)*fs[<i>,e]*ts[e,<k>]               # 8    
    HSS[ia,kc] += -ints[<i>,m,<k>,e]*ts[e,m]*(a==c)         # 5 
    HSS[ia,kc] += -0.5*(a==c)*ints[<i>,m,e,f]*td[e,f,<k>,m] # 11 
    HSS[ia,kc] += -(a==c)*ints[<i>,m,e,f]*ts[e,m]*ts[f,<k>] # 14 
    
From reference [2] table I

W<sub>ajib</sub> = g<sub>ajib</sub> + t<sup>e</sup><sub>i</sub>g<sub>ajeb</sub> - t<sup>a</sup><sub>m</sub>g<sub>mjib</sub> - (t<sup>ea</sup><sub>im</sub> + t<sup>e</sup><sub>i</sub>t<sup>a</sup><sub>m</sub>)g<sub>mjeb</sub>

We can re-index (j->k, b->c) as W<sub>**akic**</sub> = g<sub>**akic**</sub> + t<sup>e</sup><sub>**i**</sub>g<sub>**ak**e**c**</sub> - t<sup>**a**</sup><sub>m</sub>g<sub>m**ki**c</sub> - (t<sup>e**a**</sup><sub>**i**m</sub> + t<sup>e</sup><sub>**i**</sub>t<sup>**a**</sup><sub>m</sub>)g<sub>m**k**e**c**</sub>

    HSS[ia,kc] += ints[<a>,<k>,<i>,<c>]                     # 3
    HSS[ia,kc] += ints[<a>,<k>,e,<c>]*ts[e,<i>]             # 6
    HSS[ia,kc] += -ints[m,<k>,<i>,<c>]*ts[<a>,m]            # 7    
    HSS[ia,kc] += -ints[m,<k>,e,<c>]*td[e,<a>,<i>,m]         # 12
    HSS[ia,kc] += -ints[m,<k>,e,<c>]*ts[e,<i>]*ts[<a>,m]    # 13
    
### H<sub>SD</sub>
From reference [2] table I

F<sub>ia</sub> = f<sub>ia</sub> + t<sup>e</sup><sub>m</sub>g<sub>imae</sub>

We can re-index ([i->l, a->d][i->k,d->c]) f<sub>**ld**</sub>&delta;<sub>**ik**</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>m</sub>g<sub>**k**m**c**e</sub>&delta;<sub>**il**</sub>&delta;<sub>**ad**</sub>

     HSD[ia,kcld] += (i==k)*(a==c)*fs[<l>,<d>]                 # 16
     HSD[ia,kcld] += ints[<k>,m,<c>,e]*ts[e,m]*(i==l)*(a==d)   # 21 

From reference [2] table I
W<sub>aibc</sub> = g<sub>aibc</sub> - t<sup>a</sup><sub>m</sub>g<sub>mibc</sub>

We can re-index (i->l,b->c,c->d) W<sub>**alcd**</sub> = g<sub>**alcd**</sub>&delta;<sub>**ik**</sub> - t<sup>**a**</sup><sub>m</sub>g<sub>m**lcd**</sub>&delta;<sub>**ik**</sub>

     HSD[ia,kcld] += 0.5*ints[<a>,<l>,<c>,<d>]*(i==k)          # 17
     HSD[ia,kcld] += -0.5*ints[m,<l>,<c>,<d>]*ts[<a>,m]*(i==k) # 20  

There is disagreement between reference [2] and [Coupled-cluster calculations of nuclear magnetic resonance chemical shifts](www2.chemia.uj.edu.pl/~migda/Literatura/pdf/JCP03561.pdf) we have taken reference [3] which agrees with coding given later **\***

From reference [3] table III
W<sub>mnie</sub> = g<sub>mnie</sub> + t<sup>f</sup><sub>i</sub>g<sub>mnfe</sub>
    
We can re-index (m->k,n->l,e->d,f->e) W<sub>klid</sub> = g<sub>**klid**</sub> </sub>&delta;<sub>**ac**</sub>+ 0.5 t<sup>e</sup><sub>**i**</sub>g<sub>**kl**e**d**</sub></sub>&delta;<sub>**ac**</sub>

     HSD[ia,kcld] += -0.5*ints[<k>,<l>,<i>,<d>]*(a==c)         # 18
     HSD[ia,kcld] += -0.5*ints[<k>,<l>,e,<d>]*ts[e,<i>]*(a==c) # 19     
 
**\***
Psi4numpy has
```python
def build_Wmnie():
    """Eqn 7 from [Gauss:1995:3561], Table III (b)"""
    Wmnie = ccsd.get_MO('ooov').copy()
    Wmnie += np.einsum('if,mnfe->mnie', t1, ccsd.get_MO('oovv'))
    return Wmnie
```
This then is (f->e, e->d,m->k,n->l) t<sup>e</sup><sub>i</sub>g<sub>kled</sub> which agrees with last reference [3]. Finally pyscf has
```python
def Wooov(t1, t2, eris):
    Wmnie = einsum('if,mnfe->mnie', t1, eris.oovv)
    Wmnie += eris.ooov
    return Wmnie
```
This is the same as psi4numpy.

### H<sub>DS</sub>

From reference [2] table I 
W<sub>iajk</sub> = g<sub>iajk</sub> + *P*(ij)t<sup>ae</sup><sub>km</sub>g<sub>imje</sub> + 0.5&tau;<sup>ef</sup><sub>jk</sub>g<sub>iaef</sub> + t<sup>a</sup><sub>m</sub>W<sub>imjk</sub>  - *P*(ij) t<sup>e</sup><sub>j</sub> (g<sub>iaek</sub> - t<sup>af</sup><sub>mk</sub>g<sub>imef</sub> ) + t<sup>ae</sup><sub>jk</sub>F<sub>ie</sub>

+ g<sub>iajk</sub> re-index (i->k,j->i,k->j) g<sub>**kaij**</sub>&delta;<sub>**bc**</sub>

      HDS[iajb,kc] += (b==c)*ints[k,a,i,j]                            HDS[iajb,kc] -= (a==c)*ints[k,b,i,j]             
      
+*P*(ij)t<sup>ae</sup><sub>km</sub>g<sub>imje</sub> re-index  (n->m) *P*(ij) t<sup>**a**e</sup><sub>**j**m</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**bc**</sub> or *P*(ij) t<sup>e**a**</sup><sub>m**j**</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**bc**</sub>

      HDS[iajb,kc] += (b==c)*ints[k,m,i,e]*td[e,a,m,j] - \            HDS[iajb,kc] -= (a==c)*ints[k,m,i,e]*td[e,b,m,j] + \ 
                      (b==c)*ints[k,m,j,e]*td[e,a,m,i]                                (a==c)*ints[k,m,j,e]*td[e,b,m,i]   

+0.5&tau;<sup>ef</sup><sub>jk</sub>g<sub>iaef</sub> is  +0.5<sup>ef</sup><sub>jk</sub>g<sub>iaef</sub> + t<sup>e</sup><sub>j</sub>t<sup>f</sup><sub>k</sub>g<sub>iaef</sub> re-index as 0.5<sup>ef</sup><sub>**ij**</sub>g<sub>**ka**ef</sub>&delta;<sub>**bc**</sub> +  t<sup>e</sup><sub>**i**</sub>t<sup>f</sup><sub>**j**</sub>g<sub>**ka**ef</sub>&delta;<sub>**ac**</sub>
 
      HDS[iajb,kc] += 0.5*(b==c)*ints[k,a,e,f]*td[e,f,i,j]            HDS[iajb,kc] -= 0.5*(a==c)*ints[k,b,e,f]*td[e,f,i,j] 
      HDS[iajb,kc] += (b==c)*ts[e,i]*ts[f,j]*ints[k,a,e,f]            HDS[iajb,kc] -= (a==c)*ts[e,i]*ts[f,j]*ints[k,b,e,f]
            
-*P*(ij) t<sup>e</sup><sub>j</sub>g<sub>iaek</sub> re-index as -t<sup>e</sup><sub>**i**</sub>g<sub>**ka**e**j**</sub>&delta;<sub>**bc**</sub> + t<sup>e</sup><sub>**j**</sub>g<sub>**ka**e**i**</sub>&delta;<sub>**bc**</sub>

      HDS[iajb,kc] += (b==c)*ints[k,a,e,j]*ts[e,i]                   HDS[iajb,kc] += (b==c)*ints[k,a,e,i]*ts[e,j]
                     -(b==c)*ints[k,a,e,i]*ts[e,j]                                  -(b==c)*ints[k,a,e,j]*ts[e,i]
 
We have a difference of sign here reference \[3] gives this term as +*P* t<sup>e</sup><sub>i</sub>g<sub>mbej</sub> which re-indexed is +*P* t<sup>e</sup><sub>i</sub>g<sub>kaej</sub> whereas reference [2] has the minus sign.


+*P*(ij) t<sup>e</sup><sub>j</sub> t<sup>af</sup><sub>mk</sub>g<sub>imef</sub> re-index as t<sup>e</sup><sub>**i**</sub> t<sup>**a**f</sup><sub>m**j**</sub>g<sub>**k**mef</sub>&delta;<sub>**bc**</sub> or **-** t<sup>e</sup><sub>**i**</sub> t<sup>f**a**</sup><sub>m**j**</sub>g<sub>**k**mef</sub>&delta;<sub>**bc**</sub>

       HDS[iajb,kc] += (b==c)*ts[e,i]*td[f,a,m,j]*ints[k,m,e,f] - \  HDS[iajb,kc] += (b==c)*ts[e,j]*td[f,a,m,i]*ints[k,m,e,f] - \ 
                       (a==c)*ts[e,i]*td[f,b,m,j]*ints[k,m,e,f]                      (a==c)*ts[e,j]*td[f,b,m,i]*ints[k,m,e,f]
      
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

Reference [2] for W<sub>klij</sub> = g<sub>klij</sub> + *P*(ij) t<sup>e</sup><sub>j</sub>g<sub>klie</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub>g<sub>klef</sub> 
-t<sup>a</sup><sub>m</sub>g<sub>imjk</sub> re-index as - t<sup>**a**</sup><sub>m</sub>g<sub>**k**m**ij**</sub>&delta;<sub>**bc**</sub>
 
                     HDS[iajb,kc] += (a==c)*ints[k,m,i,j]*ts[b,m] - \
                                     (b==c)*ints[k,m,i,j]*ts[a,m] 
                                     

-t<sup>a</sup><sub>m</sub> *P*(ij) t<sup>e</sup><sub>j</sub>g<sub>klie</sub> re-indexing to W<sub>imjk</sub>, - t<sup>a</sup><sub>m</sub> *P*(ij) t<sup>e</sup><sub>k</sub>g<sub>imje</sub> and re-indexing (i->k,j->i,k->j) gives -t<sup>a</sup><sub>m</sub> *P*(ij) t<sup>e</sup><sub>j</sub>g<sub>kmie</sub> and swapping g-index

+t<sup>**a**</sup><sub>m</sub> *P*(ij) t<sup>e</sup><sub>**j**</sub>g<sub>**k**me**i**</sub> &delta;<**bc**</sub>

                     (b==c)*ts[e,j]*ts[a,m]*ints[k,m,e,i] # 39
                    -(b==c)*ts[e,i]*ts[a,m]*ints[k,m,e,j]                   

                    HDS[iajb,kc] -= (a==c)*ts[e,j]*ts[b,m]*ints[k,m,e,i]  +\
                                    (a==c)*ts[e,i]*ts[b,m]*ints[k,m,e,j] 
                                    
-0.5t<sup>**a**</sup><sub>m</sub>&tau;<sup>ef</sup><sub>ij</sub>g<sub>klef</sub> is -0.5t<sup>**a**</sup><sub>m</sub>t<sup>ef</sup><sub>ij</sub>g<sub>klef</sub> - t<sup>**a**</sup><sub>m</sub>t<sup>e</sup><sub>i</sub><sup>f</sup><sub>j</sub>g<sub>klef</sub>
re-index to W<sub>imjk</sub>, -0.5t<sup>**a**</sup><sub>m</sub>t<sup>ef</sup><sub>jk</sub>g<sub>imef</sub> then re-indexing (i->k,j->i,k->j) gives -0.5t<sup>**a**</sup><sub>m</sub>t<sup>ef</sup><sub>ij</sub>g<sub>kmef</sub>&delta;<sub>**bc**</sub>

                     -0.5*(b==c)*ts[a,m]*td[e,f,i,j]*ints[k,m,e,f]  # 43
                      0.5*(a==c)*ts[b,m]*td[e,f,i,j]*ints[k,m,e,f] 
                      
-t<sup>**a**</sup><sub>m</sub>t<sup>e</sup><sub>i</sub><sup>f</sup><sub>j</sub>g<sub>klef</sub> re-index as -t<sup>e</sup><sub>j</sub><sup>f</sup><sub>k</sub>g<sub>imef</sub> and the re-indexing (i->k,j->i,k->j) gives -t<sup>**a**</sup><sub>m</sub>t<sup>e</sup><sub>i</sub><sup>f</sup><sub>j</sub>g<sub>kmef</sub> and swapping index on g-tensor gives
+t<sup>**a**</sup><sub>m</sub>t<sup>e</sup><sub>**i**</sub><sup>f</sup><sub>**j**</sub>g<sub>m**k**ef</sub>&delta;<sub>**bc**</sub>

                      HDS[iajb,kc] += (b==c)*ts[e,i]*ts[a,m]*ts[f,j]*ints[m,k,e,f] - \
                                      (a==c)*ts[e,i]*ts[b,m]*ts[f,j]*ints[m,k,e,f]  # 51
 
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


                                      
                                      

