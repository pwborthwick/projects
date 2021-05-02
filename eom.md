# EOM-CCSD
The paper [Simplified methods for equation-of-motion coupled-cluster excited state calculations - Steven R. Gwaltney, Marcel Nooijen, Rodney J. Bartlett](https://notendur.hi.is/agust/rannsoknir/papers/cpl248-189-96.pdf) gives the following equations for the partitioning of the EOM-CCSD hamiltonian 

H<sub>SS</sub><sup>**a**</sup><sub>**i**</sub> = \[*P*(ij)F<sub>**a**c</sub> &delta;<sub>**i**k</sub> - *P*(ab) F<sub>k**i**</sub> &delta;<sub>**a**c</sub> + W<sub>**a**k**i**c</sub>] r<sup>ck</sup>

H<sub>SD</sub><sup>**a**</sup><sub>**i**</sub> = \[F<sub>ld</sub> &delta;<sub>**i**k</sub>&delta;<sub>**a**c</sub> + 0.5  W<sub>**a**ldc</sub> &delta;<sub>**i**k</sub> - 0.5 W<sub>lk**i**d</sub> &delta;<sub>**a**c</sub>] r<sup>lkcd</sup>

H<sub>DS</sub><sup>**ab**</sup><sub>**ij**</sub> = \[*P*(ab) W<sub>k**aij**</sub> &delta;<sub>**b**c</sub> + *P*(ij) W<sub>**ab**c**j**</sub> &delta;<sub>**i**k</sub> + *P*(ab) W<sub>**b**kec</sub> t<sup>**a**e</sup><sub>**ij**</sub> - *P*(ij) W<sub>mk**j**c</sub> t<sup>**ab**</sup><sub>**i**m]</sub>] r<sup>kc</sup> 

H<sub>DD</sub><sup>**ab**</sup><sub>**ij**</sub> = \[*P*(ab) F<sub>**b**c</sub> &delta;<sub>**j**k</sub> &delta;<sub>**i**l</sub> &delta;<sub>**a**d</sub> - F<sub>k**j**</sub></sub> &delta;<sub>**a**d</sub> &delta;<sub>**i**l</sub> &delta;<sub>**b**c</sub> + 0.5W<sub>**ab**cd</sub> &delta;<sub>**i**k</sub> &delta;<sub>**j**l</sub> + 0.5W<sub>kl**ij**</sub> &delta;<sub>**a**c</sub > &delta;<sub>**b**d</sub> + *P*(ij)*P*(ab) W<sub>**a**k**i**c</sub> &delta;<sub>**j**l</sub> &delta;<sub>**b**d</sub> - 0.5W<sub>lkec</sub>t<sup>e**b**</sup><sub>**ij**</sub> &delta;<sub>**a**d</sub> + 0.5W<sub>mkdc</sub>t<sup>**ab**</sup><sub>**j**m</sub> &delta;<sub>**i**l</sub>] r<sup>lkcd</sup>

*(Einstein summation implied on repeated indices)*
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

+ F<sub>a**c**</sub> = *f*<sub>a**c**</sub> &delta;<sub>i**k**</sub> - t<sub>a</sub><sup>m</sup> *f*<sub>m**c**</sub>&delta;<sub>i**k**</sub> + t<sup>em</sup> g<sub>mae**c**</sub> &delta;<sub>i**k**</sub> - 0.5 t<sup>e</sup><sub>a</sub><sup>nm</sup> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub> - t<sup>e</sup><sub>n</sub> t<sup>a</sup><sub>m</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>

    + [1]  +*f*<sub>a**c**</sub> &delta;<sub>i**k**</sub>
    + [2]  -t<sub>a</sub><sup>m</sup> *f*<sub>m**c**</sub> &delta;<sub>i**k**</sub>
    + [3]  +t<sup>em</sup>g<sub>mae**c**</sub> &delta;<sub>i**k**</sub> 
    + [4]  -0.5 t<sup>e</sup><sub>a</sub><sup>nm</sup> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub> 
    + [5]  -t<sup>e</sup><sub>n</sub> t<sup>a</sup><sub>m</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>
   
+ F<sub>**k**i</sub> = *f*<sub>**k**i</sub> &delta;<sub>a**c**</sub> + t<sup>e</sup><sub>i</sub> *f*<sub>**k**e</sub> &delta;<sub>a**c**</sub> + t<sup>e</sup><sub>m</sub> g<sub>**k**mie</sub> &delta;<sub>a**c**</sub> + 0.5 t<sup>ef</sup><sub>i</sub><sup>m</sup> g<sub>**k**mef</sub> &delta;<sub>a**c**</sub> + t<sup>e</sup><sub>i</sub> t<sup>f</sup><sub>m</sub> g<sub>**k**mef</sub> &delta;<sub>a**c**</sub> 

    + [6]  +*f*<sub>**k**i</sub> &delta;<sub>a**c**</sub>
    + [7]  +t<sup>e</sup><sub>i</sub> *f*<sub>**k**e</sub> &delta;<sub>a**c**</sub>
    + [8]  +t<sup>e</sup><sub>m</sub> g<sub>**k**mie</sub> &delta;<sub>a**c**</sub>
    + [9]  +0.5 t<sup>ef</sup><sub>i</sub><sup>m</sup> g<sub>**k**mef</sub> &delta;<sub>a**c**</sub>
    + [10] +t<sup>e</sup><sub>i</sub> t<sup>f</sup><sub>m</sub> g<sub>**k**mef</sub> &delta;<sub>a**c**</sub>

+ W<sub>a**k**i**c**</sub> = g<sub>a**k**i**c**</sub> + t<sup>e</sup><sub>i</sub> g<sub>a**k**e**c**</sub> - t<sup>a</sup><sub>m</sub> g<sub>m**k**i**c**</sub> - (t<sup>ea</sup><sub>im</sub> + t<sup>e</sup><sub>i</sub> t<sup>a</sup><sub>m</sub>) g<sub>m**k**e**c**</sub>

    + [11] +g<sub>a**k**i**c**</sub>
    + [12] +t<sup>e</sup><sub>i</sub> g<sub>a**k**e**c**</sub>
    + [13] -t<sup>a</sup><sub>m</sub>  g<sub>m**k**i**c**</sub>
    + [14] -t<sup>ea</sup><sub>im</sub>g<sub>m**k**e**c**</sub>
    + [15] -t<sup>e</sup><sub>i</sub> t<sup>a</sup><sub>m</sub> g<sub>m**k**e**c**</sub>
 
- - -
#### H<sub>SD</sub>
+ F<sub>**ld**</sub> = f<sub>**ld**</sub> &delta;<sub>i**k**</sub>&delta;<sub>**c**</sub> + t<sup>e</sup><sub>m</sub> g<sub>**k**m**c**e</sub> &delta;<sub>i**l**</sub>&delta;<sub>a**d**</sub>

    + [16] +f<sub>**ld**</sub> &delta;<sub>i**k**</sub>&delta;<sub>a**c**</sub>
    + [17] +t<sup>e</sup><sub>m</sub> g<sub>**k**m**c**e</sub> &delta;<sub>i**l**</sub>&delta;<sub>a**d**</sub>
    
+  W<sub>a**ldc**</sub> = g<sub>a**ldc**</sub>  &delta;<sub>i**k**</sub> - t<sup>a</sup><sub>m</sub>g<sub>m**ldc**</sub> &delta;<sub>i**k**</sub>

    + [18] +g<sub>a**ldc**</sub> &delta;<sub>a**k**</sub> 
    + [19] -t<sup>a</sup><sub>m</sub> g<sub>m**ldc**</sub> &delta;<sub>i**k**</sub>

+ W<sub>**lk**i**d**</sub> = g<sub>**lk**i**d**</sub> </sub>&delta;<sub>a**c**</sub> + t<sup>e</sup><sub>i</sub> g<sub>**lk**e**d**</sub></sub> &delta;<sub>a**c**</sub>

    + [20] +g<sub>**lk**i**d**</sub> </sub>&delta;<sub>a**c**</sub>
    + [21] +t<sup>e</sup><sub>i</sub> g<sub>**lk**e**d**</sub></sub> &delta;<sub>a**c**</sub>
    
*There is disagreement between reference [2] and [Coupled-cluster calculations of nuclear magnetic resonance chemical shifts](www2.chemia.uj.edu.pl/~migda/Literatura/pdf/JCP03561.pdf) we have taken reference [3] which agrees with coding in psi4numpy/pyscf. Reference [2] has g<sub>kild</sub> + t<sup>d</sup><sub>l</sub>g<sub>kied</sub> and reference [3] g<sub>lkid</sub> + t<sup>e</sup><sub>i</sub>g<sub>lked</sub>*

- - -
#### H<sub>DS</sub>
+ W<sub>**k**aij</sub> = g<sub>**k**aij</sub> &delta;<sub>b**c**</sub> + *P*(ij) t<sup>ea</sup><sub>mj</sub> g<sub>**k**mie</sub> &delta;<sub>b**c**</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub> g<sub>**k**aef</sub> &delta;<sub>b**c**</sub> - t<sup>a</sup><sub>m</sub> W<sub>**k**mij</sub> &delta;<sub>b**c**</sub> +*P*(ij) t<sup>e</sup><sub>i</sub> (g<sub>**k**aej</sub> - t<sup>af</sup><sub>mj</sub> g<sub>**k**mef</sub>) &delta;<sub>b**c**</sub> +t<sup>ea</sup><sub>ij</sub> F<sub>**k**e</sub> &delta;<sub>b**c**</sub>

    + [22] +g<sub>**k**aij</sub> &delta;<sub>b**c**</sub>
    + [23] +*P*(ij) t<sup>ea</sup><sub>mj</sub> g<sub>**k**mie</sub> &delta;<sub>*b*c**</sub>
    + [--] +0.5&tau;<sup>ef</sup><sub>ij</sub> g<sub>**k**aef</sub> &delta;<sub>b**c**</sub> 
        + [24] +0.5t<sup>ef</sup><sub>ij</sub> g<sub>**k**aef</sub> &delta;<sub>b**c**</sub> 
        + [25] +t<sup>e</sup><sub>i</sub> t<sup>f</sup><sub>j</sub> g<sub>**k**aef</sub> &delta;<sub>**bc**</sub>
    + [--] -t<sup>a</sup><sub>m</sub>W<sub>**k**mij</sub> &delta;<sub>b**c**</sub>
        + W<sub>**k**mij</sub> =  g<sub>**k**mij</sub> + *P*(ij) t<sup>e</sup><sub>j</sub> g<sub>**k**mie</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub> g<sub>**k**mef</sub> 
            + [26] -t<sup>a</sup><sub>m</sub> g<sub>**k**mij</sub> &delta;<sub>b**c**</sub>
            + [27] +t<sup>a</sup><sub>m</sub> *P*(ij) t<sup>e</sup><sub>j</sub> g<sub>**k**mei</sub> &delta;<sub>b**c**</sub>
            + [--] -0.5&tau;<sup>a</sup><sub>m</sub> t<sup>ef</sup><sub>ij</sub> g<sub>**k**mef</sub> &delta;<sub>b**c**</sub>
                + [28] -0.5t<sup>a</sup><sub>m</sub> t<sup>ef</sup><sub>ij</sub> g<sub>**k**mef</sub> &delta;<sub>b**c**</sub>
                + [29] -t<sup>a</sup><sub>m</sub> t<sup>e</sup><sub>i</sub> t<sup>f</sup><sub>j</sub> g<sub>**k**mef</sub> &delta;<sub>b**c**</sub>
    + [--] +*P*(ij) t<sup>e</sup><sub>i</sub> (g<sub>**k**aej</sub> - t<sup>af</sup><sub>mj</sub> g<sub>**k**mef</sub>) &delta;<sub>b**c**</sub>
        + [30] +*P*(ij) t<sup>e</sup><sub>i</sub> g<sub>**k**aej</sub> &delta;<sub>b**c**</sub> 
        + [31] -*P*(ij) t<sup>e</sup><sub>i</sub> t<sup>af</sup><sub>mj</sub>g<sub>**k**mef</sub> &delta;<sub>b**c**</sub> 
    + [--] +t<sup>ea</sup><sub>ij</sub> F<sub>**k**e</sub> &delta;<sub>b**c**</sub> = t<sup>ea</sup><sub>ij</sub> *f*<sub>**k**e</sub> &delta;<sub>**bc**</sub> + t<sup>ea</sup><sub>ij</sub> t<sup>f</sup><sub>m</sub> g<sub>**k**mef</sub> &delta;<sub>b**c**</sub>
        + [32] +t<sup>ea</sup><sub>ij</sub> *f*<sub>**k**e</sub> &delta;<sub>b**c**</sub>
        + [48] +t<sup>ea</sup><sub>ij</sub> t<sup>f</sup><sub>m</sub> g<sub>**k**mef</sub> &delta;<sub>b**c**</sub>
                                     
+ W<sub>ab**c**j</sub> = g<sub>ab**c**j</sub> &delta;<sub>i**k**</sub>  - *P*(ab) g<sub>mb**c**f</sub> t<sup>af</sup><sub>mj</sub> &delta;<sub>i**k**</sub>  + 0.5&tau;<sup>ab</sup><sub>mn</sub> g<sub>mn**c**j</sub> &delta;<sub>i**k**</sub> + t<sup>e</sup><sub>j</sub> W<sub>ab**c**e</sub> &delta;<sub>*i*k**</sub>  - *P*(ab) t<sup>a</sup><sub>m</sub> (g<sub>mb**c**j</sub> &delta;<sub>i**k**</sub>  - t<sup>be</sup><sub>nj</sub> g<sub>mn**c**e</sub>) &delta;<sub>i**k**</sub>  - F<sub>m**c**</sub> t<sup>ab</sup><sub>mj</sub> &delta;<sub>i**k**</sub> 
 
    + [33] +g<sub>ab**c**j</sub> &delta;<sub>i**k**</sub> 
    + [34] -*P*(ab) g<sub>bm**c**e</sub> t<sup>ae</sup><sub>mi</sub> &delta;<sub>j**k**</sub>     
    + [--] +0.5&tau;<sup>ab</sup><sub>mn</sub> g<sub>mn**c**j</sub> &delta;<sub>i**k**</sub>
        + [35] +0.5t<sup>ab</sup><sub>mn</sub> g<sub>mn**c**j</sub> &delta;<sub>i**k**</sub>
        + [36] +t<sup>a</sup><sub>m</sub> t<sup>b</sup><sub>n</sub> g<sub>mn**c**j</sub> &delta;<sub>i**k**</sub>                             
    + [--] W<sub>ab**c**e</sub> = (g<sub>ab**c**e</sub> - *P*(ab) t<sup>b</sup><sub>m</sub> g<sub>am**c**e</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub> g<sub>mn**c**e</sub>) &delta;<sub>i**k**</sub>
        + [37] +t<sup>e</sup><sub>j</sub> g<sub>ab**c**e</sub> &delta;<sub>i**k**</sub>
        + [38] *P*(ab) t<sup>e</sup><sub>j</sub> t<sup>b</sup><sub>m</sub> g<sub>ma**c**e</sub> &delta;<sub>i**k**</sub>
        + [--] 0.5&tau;<sup>ab</sup><sub>mn</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>
            + [39] 0.5t<sup>e</sup><sub>j</sub> t<sup>ab</sup><sub>mn</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>
            + [40] t<sup>e</sup><sub>j</sub> t<sup>a</sup><sub>m</sub> t<sup>b</sup><sub>n</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>
    + [41] -*P*(ab) t<sup>a</sup><sub>m</sub> g<sub>mb**c**j</sub> &delta;<sub>i**k**</sub>
    + [42] -*P*(ab) t<sup>a</sup><sub>m</sub> t<sup>eb</sup><sub>nj</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>
    + [--] -F<sub>m**c**</sub> t<sup>ab</sup><sub>mj</sub> &delta;<sub>i**k**</sub> = -t<sup>ab</sup><sub>mj</sub> *f*<sub>m**c**</sub> &delta;<sub>i**k**</sub> - t<sup>ab</sup><sub>mj</sub> t<sup>e</sup><sub>n</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>
        + [47] -t<sup>ab</sup><sub>mj</sub> *f*<sub>mc</sub> &delta;<sub>i**k**</sub>
        + [49] -t<sup>ab</sup><sub>mj</sub> t<sup>e</sup><sub>n</sub> g<sub>mn**c**e</sub> &delta;<sub>i**k**</sub>
   
+ *P*(ab) W<sub>bke**c**</sub> t<sup>ae</sup><sub>ij</sub> = t<sup>ae</sup><sub>ij</sub>(g<sub>bke**c**</sub> - t<sup>b</sup><sub>m</sub> g<sub>m**k**e**c**</sub>)

    + [43] +t<sup>ae</sup><sub>ij</sub> g<sub>bke**c**</sub>
    + [44] -t<sup>ae</sup><sub>ij</sub> t<sup>b</sup><sub>m</sub> g<sub>m**k**e**c**</sub>
    
+ -*P*(ij) W<sub>m**k**j**c**</sub> t<sup>ab</sup><sub>im</sub> =  -t<sup>ab</sup><sub>im</sub>(g<sub>m**k**j**c**</sub> - t<sup>e</sup><sub>j</sub> g<sub>m**k**e**c**</sub>)
    + [45] -t<sup>ab</sup><sub>im</sub> g<sub>m**k**j**c**</sub> 
    + [46] -t<sup>ab</sup><sub>im</sub> t<sup>e</sup><sub>j</sub> g<sub>m**k**e**c**</sub>                    

- - -
#### H<sub>DD</sub>
+ +*P*(ab) F<sub>b**c**</sub> = (*f*<sub>b**c**</sub> - *f*<sub>m**c**</sub> t<sup>b</sup><sub>m</sub> + t<sup>e</sup><sub>m</sub> g<sub>mae**c**</sub> - 0.5 &tau;<sup>eb</sup><sub>mn</sub>g<sub>mne**c**</sub>) &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>&delta;<sub>a**d**</sub>
    + [50] +*f*<sub>b**c**</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>&delta;<sub>a**d**</sub>
    + [51] -*f*<sub>m**c**</sub> t<sup>b</sup><sub>m</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>&delta;<sub>a**d**</sub>
    + [52] +t<sup>e</sup><sub>m</sub>g<sub>mbe**c**</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>&delta;<sub>a**d**</sub>
    + [53] -0.5t<sup>eb</sup><sub>mn</sub> g<sub>mne**c**</sub>&delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>&delta;<sub>a**d**</sub>
    + [54] -0.5t<sup>e</sup><sub>m</sub> t<sup>a</sup><sub>n</sub> g<sub>mne**c**</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>&delta;<sub>a**d**</sub>
    
+ -*P*(ij) F<sub>**k**j</sub> = (*f*<sub>**k**j</sub> + *f*<sub>**k**e</sub> t<sup>e</sup><sub>j</sub> + t<sup>e</sup><sub>m</sub> g<sub>**k**mje</sub> + 0.5&tau;<sup>ef</sup><sub>jm</sub> g<sub>**k**mef</sub> ) &delta;<sub>a**d**</sub>&delta;<sub>i**l**</sub>&delta;<sub>bc</sub>
    + [55] +*f*<sub>**k**j</sub></sub> &delta;<sub>a**d**</sub>&delta;<sub>i**l**</sub>&delta;<sub>bc</sub>
    + [56] +*f*<sub>**k**e</sub> t<sup>e</sup><sub>j</sub> &delta;<sub>a**d**</sub>&delta;<sub>i**l**</sub>&delta;<sub>bc</sub>
    + [57] +t<sup>e</sup><sub>m</sub> g<sub>**k**mie</sub>  &delta;<sub>a**d**</sub>&delta;<sub>i**l**</sub>&delta;<sub>bc</sub>
    + [58] +0.5t<sup>ef</sup><sub>im</sub> g<sub>**k**mef</sub> &delta;<sub>a**d**</sub>&delta;<sub>i**l**</sub>&delta;<sub>bc</sub>
    + [59] +t<sup>e</sup><sub>i</sub> t<sup>f</sup><sub>m</sub> g<sub>**k**mef</sub> &delta;<sub>a**d**</sub>&delta;<sub>i**l**</sub>&delta;<sub>bc</sub>
    
+ +0.5W<sub>ab**cd**</sub> = (g<sub>ab**cd**</sub> - *P*(ab) t<sup>b</sup><sub>m</sub> g<sub>am**cd**</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub> g<sub>mn**cd**</sub>) &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>
    + [60] +0.5g<sub>ab**cd**</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>
    + [61] -0.5t<sup>b</sup><sub>m</sub> g<sub>am**cd**</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>
    + [62] +0.5t<sup>ab</sup><sub>mn</sub> g<sub>mn**cd**</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>
    + [63] t<sup>a</sup><sub>m</sub> t<sup>b</sup><sub>n</sub> g<sub>mn**cd**</sub> &delta;<sub>i**k**</sub>&delta;<sub>j**l**</sub>
    
+ +0.5W<sub>**kl**ij</sub> = (g<sub>**kl**ij</sub> + *P*(ij) t<sup>e</sup><sub>j</sub> g<sub>**kl**ie</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub> g<sub>**kl**ef</sub>) &delta;<sub>a**c**</sub>&delta;<sub>b**d**</sub>
    + [64] +0.5g<sub>**kl**ij</sub> &delta;<sub>a**c**</sub>&delta;<sub>b**d**</sub>
    + [65] +t<sup>e</sup><sub>j</sub> g<sub>**kl**ie</sub> &delta;<sub>a**c**</sub>&delta;<sub>b**d**</sub>
    + [66] +0.5t<sup>ef</sup><sub>ij</sub> g<sub>**kl**ef</sub> &delta;<sub>a**c**</sub>&delta;<sub>b**d**</sub>
    + [67] +t<sup>e</sup><sub>i</sub> t<sup>f</sup><sub>j</sub> g<sub>**kl**ef</sub> &delta;<sub>a**c**</sub>&delta;<sub>b**d**</sub>
    
+ +*P*(ij)*P*(ab) W<sub>a**k**i**c**</sub> = (g<sub>a**k**i**c**</sub> + t<sup>e</sup><sub>i</sub> g<sub>a**k**e**c**</sub> - t<sup>a</sup><sub>m</sub> g<sub>m**k**i**c**</sub> - (t<sup>ea</sup><sub>im</sub> + t<sup>e</sup><sub>i</sub> t<sup>a</sup><sub>m</sub>) g<sub>m**k**e**c<**/sub>) &delta;<sub>j**l**</sub>&delta;<sub>**d**b</sub>
    + [68] +g<sub>a**k**i**c**</sub>  &delta;<sub>j**l**</sub>&delta;<sub>**d**b</sub>
    + [69] +t<sup>e</sup><sub>i</sub> g<sub>a**k**e**c**</sub> &delta;<sub>j**l**</sub>&delta;<sub>**d**b</sub>
    + [70] -t<sup>a</sup><sub>m</sub>g<sub>m**k**i**c**</sub> &delta;<sub>j**l**</sub>&delta;<sub>**d**b</sub>
    + [71] -t<sup>ea</sup><sub>im</sub> g<sub>m**k**e**c**</sub> &delta;<sub>j**l**</sub>&delta;<sub>**d**b</sub>
    + [72] -t<sup>e</sup><sub>i</sub> t<sup>a</sup><sub>m</sub> g<sub>m**k**e**c**</sub> &delta;<sub>j**l**</sub>&delta;<sub>**d**b</sub>
  
+ -0.5W<sub>**lk**e**c**</sub> t<sup>eb</sup><sub>ij</sub> = -0.5t<sup>eb</sup><sub>ij</sub> g<sub>**lk**e**c**</sub> &delta;<sub>b**d**</sub>
    + [73] -0.5t<sup>eb</sup><sub>ij</sub>g<sub>**lk**e**c**</sub> &delta;<sub>b**d**</sub>

+ +0.5W<sub>m**kdc**</sub> t<sup>ab</sup><sub>jm</sub> = +0.5W<sub>m**kdc**</sub >t<sup>ab</sup><sub>jm</sub> &delta;<sub>i**l**</sub>
    + [74] -0.5W<sub>**k**m**cd**</sub> t<sup>ab</sup><sub>mj</sub> &delta;<sub>i**l**</sub>
