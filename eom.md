# EOM-CCSD
The paper [Simplified methods for equation-of-motion coupled-cluster excited state calculations - Steven R. Gwaltney, Marcel Nooijen, Rodney J. Bartlett](https://notendur.hi.is/agust/rannsoknir/papers/cpl248-189-96.pdf) gives the following equations for the partitioning of the EOM-CCSD hamiltonian 

H<sub>SS</sub><sup>**a**</sup><sub>**i**</sub> = *P*(ij)F<sub>**ac**</sub>r<sup>**c**</sup><sub>i</sub> &delta;<sub>**ik**</sub> - *P*(ab) F<sub>**ki**</sub>r<sup>**k**</sup><sub>**a**</sub> &delta;<sub>ac</sub> + W<sub>akic</sub>r<sup>**ck**</sup>

H<sub>SD</sub><sup>**a**</sup><sub>**i**</sub> = F<sub>ld</sub>r<sup>ld</sup><sub>**ai**</sub>&delta;<sub>**ik**</sub>&delta;<sub>**ac**</sub> + 0.5  W<sub>**a**ld**c**</sub>r<sup>ld**c**</sup><sub>**i**</sub> &delta;<sub>**ik**</sub> - 0.5 W<sub>l**ki**d</sub>r<sup>l**k**d</sup><sub>**a**</sub>&delta;<sub>**ac**</sub>

H<sub>DS</sub><sup>**ab**</sup><sub>**ij**</sub> = *P*(ab) W<sub>**kai**j</sub>r<sup>k</sup><sub>**b**</sub>&delta;<sub>**b**c</sub> + *P*(ij) W<sub>**ab**c**j**</sub>r<sup>**c**</sup><sub>**i**</sub>&delta;<sub>**ik**</sub> + *P*(ab) W<sub>**b**kec</sub> t<sup>**a**e</sup><sub>**ij**</sub>r<sup>ke**c**</sup> - *P*(ij) W<sub>mk**j**c</sub> t<sup>**ab**</sup><sub>**i**m</sub>r<sup>mkc</sup> 

H<sub>DD</sub><sup>**ab**</sup><sub>**ij**</sub> = *P*(ab) F<sub>**b**c</sub>&delta;<sub>**j**k</sub>&delta;<sub>**i**l</sub>&delta;<sub>**a**d</sub>r<sup>c</sup><sub>**aij**</sub> - F<sub>k**j**</sub></sub>&delta;<sub>**a**d</sub>&delta;<sub>**i**l</sub>&delta;<sub>**b**c</sub>r<sup>k</sup><sub>**iab**</sub> + 0.5W<sub>**ab**cd</sub>&delta;<sub>**i**k</sub>&delta;<sub>**j**l</sub>r<sup>cd</sup><sub>**ij**</sub> + 0.5W<sub>kl**ij**</sub>&delta;<sub>**a**c</sub>&delta;<sub>**b**d</sub>r<sup>**ij**</sup><sub>**ab**</sub> + *P*(ij)*P*(ab) W<sub>**a**k**i**c</sub>&delta;<sub>**j**l</sub>&delta;<sub>**b**d</sub>r<sup>kc</sup><sub>**bj**</sub> - 0.5W<sub>lke**c**</sub>t<sup>e**b**</sup><sub>**ij**</sub>&delta;<sub>**a**d</sub>r<sup>**i**kc</sup><sub>**a**</sub> + 0.5W<sub>mkd**c**</sub>t<sup>**ab**</sup><sub>**j**m</sub>&delta;<sub>**i**l</sub>r<sup>kd**c**</sup><sub>**i**</sub>

(Einstein summation implied on repeated indices)
- - -
Note g<sub>abcd</sub> = <ab||cd> = -<ba||cd> = -<ab||dc> = <ba||dc> \
&tau;<sup>ab</sup><sub>ij</sub> = t<sup>ab</sup><sub>ij</sub> + t<sup>a</sup><sub>i</sub><sup>b</sup><sub>j</sub> \
&tau;<sup>ab</sup><sub>ij</sub> = -&tau;<sup>ba</sup><sub>ij</sub> = -&tau;<sup>ab</sup><sub>ji</sub> \
P(ij) = f(ij) - f(ji)
- - - -
#### Intermediates
F<sub>me</sub> = F<sub>ov</sub> = *f*<sub>me</sub> + t<sup>f</sup><sub>n</sub><mn||ef> \
F<sub>mi</sub> = F<sub>oo</sub> = *f*<sub>mi</sub> + t<sup>e</sup><sub>i</sub>*f*<sub>me</sub> + t<sup>e</sup><sub>n</sub><mn||ie> + &tau;<sup>ef</sup><sub>in</sub><mn||ef> \
F<sub>ae</sub> = F<sub>vv</sub> = *f*<sub>ae</sub> - t<sup>e</sup><sub>a</sub>*f*<sub>me</sub> + t<sup>f</sup><sub>m</sub><am||ef> - &tau;<sup>fa</sup><sub>ma</sub><mn||fe> \
\
W<sub>mnij</sub> = W<sub>oooo</sub> = <mn||ij> + P(ij) t<sup>e</sup><sub>j</sub><mn||ie> + 0.5&tau;<sup>ef</sup><sub>ij</sub><mn||ef> \
W<sub>abef</sub> = W<sub>vvvv</sub> = <ab||ef> - P(ab) t<sup>b</sup><sub>m</sub><am||ef> + 0.5&tau;<sup>ab</sup><sub>mn</sub><mn||ef> \
W<sub>amef</sub> = W<sub>vovv</sub> = <am||ef> - t<sup>a</sup><sub>n</sub><nm||ef> \
W<sub>mnie</sub> = W<sub>ooov</sub> = <mn||ie> + t<sup>f</sup><sub>i</sub><mn||fe> \
W<sub>mbej</sub> = W<sub>ovvo</sub> = <mb||ej> + t<sup>f</sup><sub>j</sub><mb||ef> - t<sup>b</sup><sub>n</sub><mn||ej> - (t<sup>fb</sup><sub>jn</sub> + t<sup>f</sup><sub>j</sub>t<sup>b</sup><sub>n</sub>)<nm||fe> \
W<sub>mbje</sub> = W<sub>ovov</sub> = <mb||je> + t<sup>f</sup><sub>j</sub><bm||ef> - t<sup>b</sup><sub>n</sub><mn||je> - (t<sup>fb</sup><sub>jn</sub> + t<sup>f</sup><sub>j</sub>t<sup>b</sup><sub>n</sub>)<nm||ef> \
W<sub>abei</sub> = W<sub>vvvo</sub> = <ab||ei> - F<sub>me</sub>t<sup>ab</sup><sub>mi</sub> + t<sup>f</sup><sub>i</sub>W<sub>abef</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub><mn||ei> - P(ab) t<sup>af</sup><sub>mi</sub><mb||ef> - P(ab) t<sup>a</sup><sub>m</sub> {<mb||ei> - t<sup>bf</sup><sub>ni</sub><mn||ef>} \
W<sub>mbij</sub> = W<sub>ovoo</sub> = <mb||ij> - F<sub>me</sub>t<sup>be</sup><sub>ij</sub> - t<sup>b</sup><sub>n</sub>W<sub>mnij</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub><mb||ef> + P(ij) t<sup>be</sup><sub>jn</sub><mn||ie> + P(ij) t<sup>e</sup><sub>i</sub> {<mb||ej> - t<sup>bf</sup><sub>nj</sub><mn||ef>} 
- - - -
#### H<sub>SS</sub>
Equations for terms taken from  [J. Chem. Phys. 98, 7029 (1993); https://doi.org/10.1063/1.46474698, 7029© 1993 American Institute of Physics.The equation of motion coupled-clustermethod. A systematic biorthogonal approach to molecular excitation energies, transition probabilities, and excited state properties](https://www.theochem.ru.nl/files/local/jcp-98-7029-1993.pdf) 

+ F<sub>**ac**</sub> = *f*<sub>**ac**</sub> &delta;<sub>**ik**</sub> - t<sub>**a**</sub><sup>m</sup> *f*<sub>m**c**</sub>&delta;<sub>**ik**</sub> + t<sup>em</sup>g<sub>m**a**e**c**</sub>&delta;<sub>**ik**</sub> - 0.5 t<sup>e</sup><sub>**a**</sub><sup>nm</sup>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub> - t<sup>e</sup><sub>n</sub>t<sup>**a**</sup><sub>m</sub>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub>

    + [1]  +*f*<sub>**ac**</sub> &delta;<sub>**ik**</sub>
    + [2]  -t<sub>**a**</sub><sup>m</sup> *f*<sub>m**c**</sub>&delta;<sub>**ik**</sub>
    + [3]  +t<sup>em</sup>g<sub>m**a**e**c**</sub>&delta;<sub>**ik**</sub> 
    + [4]  -0.5 t<sup>e</sup><sub>**a**</sub><sup>nm</sup>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub> 
    + [5]  -t<sup>e</sup><sub>n</sub>t<sup>**a**</sup><sub>m</sub>g<sub>mn**c**e</sub>&delta;<sub>**ik**</sub>
   
+ F<sub>**ki**</sub> = *f*<sub>**ki**</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>**i**</sub> *f*<sub>**k**e</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>m</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**ac**</sub> + 0.5 t<sup>ef</sup><sub>**i**</sub><sup>m</sup>g<sub>**k**mef</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>**i**</sub>t<sup>f</sup><sub>m</sub>g<sub>**k**mef</sub>&delta;<sub>**ac**</sub> 

    + [6]  +*f*<sub>**ki**</sub>&delta;<sub>**ac**</sub>
    + [7]  +t<sup>e</sup><sub>**i**</sub> *f*<sub>**k**e</sub>&delta;<sub>**ac**</sub>
    + [8]  +t<sup>e</sup><sub>m</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**ac**</sub> 
    + [9]  +0.5 t<sup>ef</sup><sub>**i**</sub><sup>m</sup>g<sub>**k**mef</sub>&delta;<sub>**ac**</sub>
    + [10] +t<sup>e</sup><sub>**i**</sub>t<sup>f</sup><sub>m</sub>g<sub>**k**mef</sub>&delta;<sub>**ac**</sub> 

+ W<sub>**akic**</sub> = g<sub>**akic**</sub> + t<sup>e</sup><sub>**i**</sub>g<sub>**ak**e**c**</sub> - t<sup>**a**</sup><sub>m</sub>g<sub>m**ki**c</sub> - (t<sup>e**a**</sup><sub>**i**m</sub> + t<sup>e</sup><sub>**i**</sub>t<sup>**a**</sup><sub>m</sub>)g<sub>m**k**e**c**</sub>

    + [11] +g<sub>**akic**</sub>
    + [12] +t<sup>e</sup><sub>**i**</sub>g<sub>**ak**e**c**</sub>
    + [13] -t<sup>**a**</sup><sub>m</sub>g<sub>m**ki**c</sub>
    + [14] -t<sup>e**a**</sup><sub>**i**m</sub>g<sub>m**k**e**c**</sub>
    + [15] -t<sup>e</sup><sub>**i**</sub>t<sup>**a**</sup><sub>m</sub>g<sub>m**k**e**c**</sub>
 
- - -
#### H<sub>SD</sub>
+ F<sub>**ld**</sub> = f<sub>**ld**</sub>&delta;<sub>**ik**</sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>m</sub>g<sub>**k**m**c**e</sub>&delta;<sub>**il**</sub>&delta;<sub>**ad**</sub>

    + [16] +f<sub>**ld**</sub>&delta;<sub>**ik**</sub>&delta;<sub>**ac**</sub>
    + [17] +t<sup>e</sup><sub>m</sub>g<sub>**k**m**c**e</sub>&delta;<sub>**il**</sub>&delta;<sub>**ad**</sub>
    
+  W<sub>**aldc**</sub> = g<sub>**aldc**</sub>&delta;<sub>**ik**</sub> - t<sup>**a**</sup><sub>m</sub>g<sub>m**ldc**</sub>&delta;<sub>**ik**</sub>

    + [18] +g<sub>**aldc**</sub>&delta;<sub>**ik**</sub> 
    + [19] -t<sup>**a**</sup><sub>m</sub>g<sub>m**ldc**</sub>&delta;<sub>**ik**</sub>

+ W<sub>**lkid**</sub> = g<sub>**lkid**</sub> </sub>&delta;<sub>**ac**</sub> + t<sup>e</sup><sub>**i**</sub>g<sub>**lk**e**d**</sub></sub>&delta;<sub>**ac**</sub>

    + [20] +g<sub>**lkid**</sub> </sub>&delta;<sub>**ac**</sub>
    + [21] +t<sup>e</sup><sub>**i**</sub>g<sub>**lk**e**d**</sub></sub>&delta;<sub>**ac**</sub>
    
*There is disagreement between reference [2] and [Coupled-cluster calculations of nuclear magnetic resonance chemical shifts](www2.chemia.uj.edu.pl/~migda/Literatura/pdf/JCP03561.pdf) we have taken reference [3] which agrees with coding in psi4numpy/pyscf. Reference [2] has g<sub>kild</sub> + t<sup>d</sup><sub>l</sub>g<sub>kied</sub> and reference [3] g<sub>lkid</sub> + t<sup>e</sup><sub>i</sub>g<sub>lked</sub>*

- - -
#### H<sub>DS</sub>
+ W<sub>kaij</sub> = g<sub>**kaij**</sub>&delta;<sub>**bc**</sub> + *P*(ij) t<sup>e**a**</sup><sub>m**j**</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**bc**</sub> + 0.5&tau;<sup>ef</sup><sub>**ij**</sub>g<sub>**ka**ef</sub>&delta;<sub>**bc**</sub> - t<sup>**a**</sup><sub>m</sub>W<sub>**k**m**ij**</sub>&delta;<sub>**bc**</sub> + +*P*(ij) t<sup>e</sup><sub>**i**</sub> (g<sub>**ka**e**j**</sub> - t<sup>**a**f</sup><sub>m**j**</sub>g<sub>**k**mef</sub>)&delta;<sub>**bc**</sub> +t<sup>ea</sup><sub>ij</sub>F<sub>ke</sub>&delta;<sub>**bc**</sub>

    + [22] +g<sub>**kaij**</sub>&delta;<sub>**bc**</sub>
    + [23] +*P*(ij) t<sup>e**a**</sup><sub>m**j**</sub>g<sub>**k**m**i**e</sub>&delta;<sub>**bc**</sub>
    + [--] +0.5&tau;<sup>ef</sup><sub>**ij**</sub>g<sub>**ka**ef</sub>&delta;<sub>**bc**</sub> 
        + [24] +0.5t<sup>ef</sup><sub>**ij**</sub>g<sub>**ka**ef</sub>&delta;<sub>**bc**</sub> 
        + [25] +t<sup>e</sup><sub>**i**</sub>t<sup>f</sup><sub>**j**</sub>g<sub>**ka**ef</sub>&delta;<sub>**ac**</sub>
    + [--] -t<sup>a</sup><sub>m</sub>W<sub>kmij</sub>&delta;<sub>**bc**</sub>
        + W<sub>**k**m**ij**</sub> =  g<sub>**k**m**ij**</sub> + *P*(ij) t<sup>e</sup><sub>**j**</sub>g<sub>**k**m**i**e</sub> + 0.5&tau;<sup>ef</sup><sub>**ij**</sub>g<sub>****kmef</sub> 
            + [26] -t<sup>**a**</sup><sub>m</sub>g<sub>**k**m**ij**</sub> &delta;<sub>**bc**</sub>
            + [27] +t<sup>**a**</sup><sub>m</sub> *P*(ij) t<sup>e</sup><sub>**j**</sub>g<sub>**k**me**i**</sub> &delta;<sub>**bc**</sub>
            + [--] -0.5&tau;<sup>**a**</sup><sub>m</sub>t<sup>ef</sup><sub>ij</sub>g<sub>kmef</sub> &delta;<sub>**bc**</sub>
                + [28] -0.5t<sup>**a**</sup><sub>m</sub>t<sup>ef</sup><sub>ij</sub>g<sub>kmef</sub> &delta;<sub>**bc**</sub>
                + [29] -t<sup>**a**</sup><sub>m</sub>t<sup>e</sup><sub>**i**</sub>t<sup>f</sup><sub>**j**</sub>g<sub>**k**mef</sub> &delta;<sub>**bc**</sub>
    + [--] +*P*(ij) t<sup>e</sup><sub>**i**</sub> (g<sub>**ka**e**j**</sub> - t<sup>**a**f</sup><sub>m**j**</sub>g<sub>**k**mef</sub>)&delta;<sub>**bc**</sub>
        + [30] +*P*(ij) t<sup>e</sup><sub>**i**</sub>g<sub>**ka**e**j**</sub>&delta;<sub>**bc**</sub> 
        + [31] -*P*(ij) t<sup>e</sup><sub>**i**</sub> t<sup>**a**f</sup><sub>m**j**</sub>g<sub>**k**mef</sub>&delta;<sub>**bc**</sub> 
    + [--] +t<sup>ea</sup><sub>ij</sub>F<sub>ke</sub>&delta;<sub>**bc**</sub> = t<sup>ea</sup><sub>ij</sub>*f*<sub>ke</sub>&delta;<sub>**bc**</sub> + t<sup>ea</sup><sub>ij</sub>t<sup>f</sup><sub>m</sub>g<sub>kmef</sub>&delta;<sub>**bc**</sub>
        + [32] +t<sup>e**a**</sup><sub>**ij**</sub>*f*<sub>**k**e</sub>&delta;<sub>**bc**</sub>
        + [48] +t<sup>e**a**</sup><sub>**ij**</sub>t<sup>f</sup><sub>m</sub>g<sub>**k**mef</sub>&delta;<sub>**bc**</sub>
                                     
+ W<sub>abcj</sub> = g<sub>**abcj**</sub>&delta;<sub>**ik**</sub>  - *P*(ab) g<sub>m**bc**f</sub>t<sup>**a**f</sup><sub>m**j**</sub>&delta;<sub>**ik**</sub>  + 0.5&tau;<sup>**ab**</sup><sub>mn</sub>g<sub>mn**cj**</sub> &delta;<sub>**ik**</sub> + t<sup>e</sup><sub>**j**</sub>W<sub>**abc**e</sub>&delta;<sub>**ik**</sub>  - *P*(ab) t<sup>**a**</sup><sub>m</sub>(g<sub>m**bcj**</sub>&delta;<sub>**ik**</sub>  - t<sup>**b**e</sup><sub>n**j**</sub>g<sub>mn**c**e</sub>)&delta;<sub>**ik**</sub>  - F<sub>m**c**</sub>t<sup>**ab**</sup><sub>m**j**</sub>&delta;<sub>**ik**</sub> 
 
    + [33] +g<sub>**abcj**</sub>&delta;<sub>**ik**</sub> 
    + [34] -*P*(ab)g<sub>**b**m**c**e</sub>t<sup>**a**e</sup><sub>m**i**</sub>&delta;<sub>**jk**</sub>     
    + [--] +0.5&tau;<sup>**ab**</sup><sub>mn</sub>g<sub>mn**cj**</sub>&delta;<sub>**ik**</sub>
        + [35] +0.5t<sup>**ab**</sup><sub>mn</sub>g<sub>mn**cj**</sub>&delta;<sub>**ik**</sub>
        + [36] +t<sup>a</sup><sub>m</sub>t<sup>**b**</sup><sub>n</sub>g<sub>mn**cj**</sub>&delta;<sub>**ik**</sub>                             
    + [--] W<sub>**abc**e</sub> = g<sub>**abc**e</sub> - *P*(ab) t<sup>**b**</sup><sub>m</sub>g<sub>**a**m**c**e</sub> + 0.5&tau;<sup>**ab**</sup><sub>mn</sub>g<sub>mn**c**e</sub>
        + [37] +t<sup>e</sup><sub>**j**</sub>g<sub>**abc**e</sub>&delta;<sub>**ik**</sub>
        + [38] *P*(ab) t<sup>e</sup><sub>**j**</sub>t<sup>**b**</sup><sub>m</sub>g<sub>m**ac**e</sub>&delta;<sub>**ik**</sub>
        + [--] 0.5&tau;<sup>**ab**</sup><sub>mn</sub>g<sub>mn**c**e</sub> &delta;<sub>**ik**</sub>
            + [39] 0.5t<sup>e</sup><sub>**j**</sub>t<sup>**ab**</sup><sub>mn</sub>g<sub>mn**c**e</sub> &delta;<sub>**ik**</sub>
            + [40] t<sup>e</sup><sub>**j**</sub>t<sup>**a**</sup><sub>m</sub>t<sup>**b**</sup><sub>n</sub>g<sub>mn**c**e</sub> &delta;<sub>**ik**</sub>
    + [41] -*P*(ab) t<sup>**a**</sup><sub>m</sub>g<sub>m**bc**j</sub> &delta;<sub>**ik**</sub>
    + [42] -*P*(ab) t<sup>**a**</sup><sub>m</sub>t<sup>e**b**</sup><sub>n**j**</sub>g<sub>mn**c**e</sub> &delta;<sub>**ik**</sub>
    + [--] -F<sub>m**c**</sub>t<sup>**ab**</sup><sub>mj</sub> &delta;<sub>**ik**</sub> = -t<sup>**ab**</sup><sub>mj</sub>*f*<sub>m**c**</sub>&delta;<sub>**ik**</sub> - t<sup>**ab**</sup><sub>mj</sub>t<sup>e</sup><sub>n</sub>g<sub>mn**c**e</sub> &delta;<sub>**ik**</sub>
        + [47] -t<sup>**ab**</sup><sub>mj</sub>*f*<sub>mc</sub>&delta;<sub>**ik**</sub>
        + [49] -t<sup>**ab**</sup><sub>mj</sub>t<sup>e</sup><sub>n</sub>g<sub>mn**c**e</sub> &delta;<sub>**ik**</sub>
   
+ *P*(ab) W<sub>**bk**e**c**</sub> t<sup>**a**e</sup><sub>**ij**</sub> = t<sup>**a**e</sup><sub>**ij**</sub>(g<sub>**bk**e**c**</sub> - t<sup>**b**</sup><sub>m</sub>g<sub>m**k**e**c**</sub>)

    + [43] +t<sup>**a**e</sup><sub>**ij**</sub>g<sub>**bk**e**c**</sub>
    + [44] -t<sup>**a**e</sup><sub>**ij**</sub>t<sup>**b**</sup><sub>m</sub>g<sub>m**k**e**c**</sub>
    
+ -*P*(ij) W<sub>m**kjc**</sub> t<sup>**ab**</sup><sub>**i**m</sub> =  -t<sup>**ab**</sup><sub>**i**m</sub>(g<sub>m**kjc**</sub> - t<sup>e</sup><sub>**j**</sub>g<sub>m**k**e**c**</sub>)
    + [45] -t<sup>**ab**</sup><sub>**i**m</sub>g<sub>m**kjc**</sub> 
    + [46] -t<sup>**ab**</sup><sub>**i**m</sub>t<sup>e</sup><sub>**j**</sub>g<sub>m**k**e**c**</sub>                    

- - -
#### H<sub>DD</sub>
+ +*P*(ab) F<sub>**bc**</sub> = *f*<sub>**bc**</sub> - *f*<sub>m**c**</sub>t<sup>**b**</sup><sub>m</sub> + t<sup>e</sup><sub>m</sub>g<sub>m**a**e**c**</sub> - 0.5 &tau;<sup>e**b**</sup><sub>mn</sub>g<sub>mne**c**</sub> 
    + [50] +*f*<sub>**bc**</sub>&delta;<sub>ik</sub>&delta;<sub>jl</sub>&delta;<sub>ad</sub>
    + [51] -*f*<sub>m**c**</sub>t<sup>**b**</sup><sub>m</sub>&delta;<sub>ik</sub>&delta;<sub>jl</sub>&delta;<sub>ad</sub>
    + [52] +t<sup>e</sup><sub>m</sub>g<sub>m**b**e**c**</sub>&delta;<sub>ik</sub>&delta;<sub>jl</sub>&delta;<sub>ad</sub>
    + [53] -0.5t<sup>e**b**</sup><sub>mn</sub>g<sub>mne**c**</sub>&delta;<sub>ik</sub>&delta;<sub>jl</sub>&delta;<sub>ad</sub>
    + [54] -0.5t<sup>**e**</sup><sub>m</sub>t<sup>a</sup><sub>n</sub>g<sub>mne**c**</sub> <sub>ik</sub>&delta;<sub>jl</sub>&delta;<sub>bd</sub>
    
+ -*P*(ij) F<sub>kj</sub> = *f*<sub>kj</sub> + *f*<sub>ke</sub>t<sup>e</sup><sub>j</sub> + t<sup>e</sup><sub>m</sub>g<sub>kmje</sub> + 0.5&tau;<sup>ef</sup><sub>jm</sub>g<sub>kmef</sub> 
    + [55] +*f*<sub>**kj**</sub></sub>&delta;<sub>ad</sub>&delta;<sub>il</sub>&delta;<sub>bc</sub>
    + [56] +*f*<sub>**k**e</sub>t<sup>e</sup><sub>j</sub>
    + [57] +t<sup>e</sup><sub>m</sub>g<sub>kmie</sub> 
    + [58] +0.5t<sup>ef</sup><sub>im</sub>g<sub>kmef</sub> 
    + [59] +t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>m</sub>g<sub>kmef</sub> 
    
+ +0.5W<sub>abcd</sub> = g<sub>abcd</sub> - *P*(ab) t<sup>b</sup><sub>m</sub>g<sub>amcd</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub>g<sub>mncd</sub>
    + [60] +0.5g<sub>abcd</sub>
    + [61] -0.5t<sup>b</sup><sub>m</sub>g<sub>amcd</sub>
    + [62] +0.5t<sup>ab</sup><sub>mn</sub>g<sub>mncd</sub>
    + [63] t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>g<sub>mncd</sub>
    
+ +0.5W<sub>klij</sub> = g<sub>klij</sub> + *P*(ij) t<sup>e</sup><sub>j</sub>g<sub>klie</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub>g<sub>klef</sub>
    + [64] +0.5g<sub>klij</sub> 
    + [65] +t<sup>e</sup><sub>j</sub>g<sub>klie</sub>
    + [66] +0.5t<sup>ef</sup><sub>ij</sub>g<sub>klef</sub>
    + [67] +t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>j</sub>g<sub>klef</sub>
    
+ +*P*(ij)*P*(ab) W<sub>akic</sub> = g<sub>akic</sub> + t<sup>e</sup><sub>i</sub>g<sub>akec</sub> - t<sup>a</sup><sub>m</sub>g<sub>mkic</sub> - (t<sup>ea</sup><sub>im</sub> + t<sup>e</sup><sub>i</sub>t<sup>a</sup><sub>m</sub>)g<sub>mkec</sub>
    + [68] +g<sub>akic</sub>
    + [69] +t<sup>e</sup><sub>i</sub>g<sub>akec</sub>
    + [70] -t<sup>a</sup><sub>m</sub>g<sub>mkic</sub>
    + [71] -t<sup>ea</sup><sub>im</sub> g<sub>mkec</sub>
    + [72] -t<sup>e</sup><sub>i</sub>t<sup>a</sup><sub>m</sub>g<sub>mkec</sub>
  
+ -0.5W<sub>lkec</sub>t<sup>eb</sup><sub>ij</sub> = -0.5t<sup>eb</sup><sub>ij</sub>g<sub>lkec</sub>
    + [73] -0.5t<sup>eb</sup><sub>ij</sub>g<sub>lkec</sub>

+ +0.5W<sub>mkdc</sub>t<sup>ab</sup><sub>jm</sub> = +0.5W<sub>mkdc</sub>t<sup>ab</sup><sub>jm</sub>
    + [74] -0.5W<sub>kmcd</sub>t<sup>ab</sup><sub>mj</sub>
