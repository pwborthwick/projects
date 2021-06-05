## CCSD Lambda Equations
These equations are taken from [Gauss and Stanton J. Chem. Phys., Vol. 103, No. 9, 1 September 1995](http://www2.chemia.uj.edu.pl/~migda/Literatura/pdf/JCP03561.pdf)\
**Lambda intermediates**
+ F<sup>&Lambda;</sup><sub>ae</sub> = F<sub>ae</sub> - 0.5t<sup>a</sup><sub>m</sub>F<sub>me</sub>
+ F<sup>&Lambda;</sup><sub>mi</sub> = F<sub>mi</sub> + 0.5t<sup>e</sup><sub>i</sub>F<sub>me</sub>
+ F<sup>&Lambda;</sup><sub>me</sub> = F<sub>mi</sub>

+ + F<sub>ae</sub> = *f*<sub>ae</sub> - t<sup>a</sup><sub>m</sub>*f*<sub>me</sub> + t<sup>f</sup><sub>m</sub> g<sub>amef</sub> - 0.5 &tau;<sup>af</sup><sub>mn</sub> g<sub>mnef</sub> 
+ + F<sub>mi</sub> = *f*<sub>mi</sub> + t<sup>e</sup><sub>i</sub>*f*<sub>me</sub> + t<sup>e</sup><sub>n</sub> g<sub>mnie</sub> + 0.5 &tau;<sup>ef</sup><sub>in</sub> g<sub>mnef</sub> 
+ + F<sub>me</sub> = *f*<sub>me</sub> + t<sup>f</sup><sub>n</sub> g<sub>mnef</sub>

+ W<sup>&Lambda;</sup><sub>mnij</sub> = W<sub>mnij</sub> + 0.25&tau;<sup>ef</sup><sub>ij</sub>g<sub>mnef</sub>
+ W<sup>&Lambda;</sup><sub>abef</sub> = W<sub>abef</sub> + 0.25&tau;<sup>ab</sup><sub>mn</sub>g<sub>mnef</sub>
+ W<sup>&Lambda;</sup><sub>mbej</sub> = W<sub>mbej</sub> - 0.5&tau;<sup>fb</sup><sub>jn</sub>g<sub>mnef</sub>

+ + W<sub>mnij</sub> = g<sub>mnij</sub> + P(ij) t<sup>e</sup><sub>j</sub> g<sub>mnie</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub> g<sub>mnef</sub> 
+ + W<sub>abef</sub> = g<sub>abef</sub> - P(ab) t<sup>b</sup><sub>m</sub> g<sub>amef</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub> g<sub>mnef</sub> 
+ + W<sub>mbej</sub> = W<sub>(ovvo)</sub> = g<sub>mbej</sub> + t<sup>f</sup><sub>j</sub> g<sub>mbef</sub> - t<sup>b</sup><sub>n</sub> g<sub>mnej</sub> - (t<sup>fb</sup><sub>jn</sub> + t<sup>f</sup><sub>j</sub>t<sup>b</sup><sub>n</sub>) g<sub>nmfe</sub> \

+ W<sup>&Lambda;</sup><sub>mnie</sub> = W<sub>mnie</sub>
+ W<sup>&Lambda;</sup><sub>amef</sub> = W<sub>amef</sub>
+ W<sup>&Lambda;</sup><sub>mbij</sub> = W<sub>mbij</sub>
+ W<sup>&Lambda;</sup><sub>abei</sub> = W<sub>abei</sub>

+ + W<sub>mnie</sub> = g<sub>mnie</sub> + t<sup>f</sup><sub>i</sub> g<sub>mnfe</sub> 
+ + W<sub>amef</sub> = g<sub>amef</sub> - t<sup>a</sup><sub>n</sub> g<sub>nmef</sub>
+ + W<sub>mbij</sub> = g<sub>mbij</sub> - F<sub>me</sub>t<sup>be</sup><sub>ij</sub> - t<sup>b</sup><sub>n</sub>W<sub>mnij</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub> g<sub>mbef</sub> + P(ij) t<sup>be</sup><sub>jn</sub> g<sub>mnie</sub> + P(ij) t<sup>e</sup><sub>i</sub> {g<sub>mbej</sub> - t<sup>bf</sup><sub>nj</sub> g<sub>mnef</sub>} 
+ + W<sub>abei</sub> = g<sub>abei</sub> - F<sub>me</sub>t<sup>ab</sup><sub>mi</sub> + t<sup>f</sup><sub>i</sub>W<sub>abef</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub> g<sub>mnei</sub> - P(ab) t<sup>af</sup><sub>mi</sub> g<sub>mbef</sub> - P(ab) t<sup>a</sup><sub>m</sub> {g<sub>mbei</sub> - t<sup>bf</sup><sub>ni</sub> g<sub>mnef</sub>} 

**Three-body Terms**
+ G<sub>ae</sub> = -0.5t<sup>ef</sup><sub>mn</sub>&lambda;<sup>af</sup><sub>mn</sub>
+ G<sub>mi</sub> = 0.5t<sup>ef</sup><sub>mn</sub>&lambda;<sup>ef</sup><sub>in</sub>

**Lambda Equations**\
**&Lambda;<sub>1</sub>** (&lambda;<sup>a</sup><sub>i</sub>) = F<sup>&Lambda;</sup><sub>ia</sub> + &lambda;<sup>e</sup><sub>i</sub>F<sup>&Lambda;</sup><sub>ea</sub> - &lambda;<sup>a</sup><sub>m</sub>F<sup>&Lambda;</sup><sub>im</sub> + &lambda;<sup>e</sup><sub>m</sub>W<sup>&Lambda;</sup><sub>ieam</sub> + 0.5&lambda;<sup>ef</sup><sub>im</sub>W<sup>&Lambda;</sup><sub>efam</sub> - 0.5&lambda;<sup>ae</sup><sub>mn</sub>W<sup>&Lambda;</sup><sub>iemn</sub> - G<sub>ef</sub>W<sup>&Lambda;</sup><sub>eifa</sub> - G<sub>mn</sub>W<sub>mina</sub>\
**&Lambda;<sub>2</sub>** (&lambda;<sup>ab</sup><sub>ij</sub>) = g<sub>ijab</sub>  - P(ij)&lambda;<sup>ab</sup><sub>im</sub>F<sup>&Lambda;</sup><sub>jm</sub> + 0.5&lambda;<sup>ab</sup><sub>mn</sub>W<sup>&Lambda;</sup><sub>ijmn</sub> + 0.5&lambda;<sup>ef</sup><sub>ij</sub>W<sup>&Lambda;</sup><sub>efab</sub> + P(ij)&lambda;<sup>e</sup><sub>i</sub>W<sup>&Lambda;</sup><sub>ejab</sub> - P(ab)&lambda;<sup>a</sup><sub>m</sub>W<sup>&Lambda;</sup><sub>ijmb</sub> + P(ij)P(ab)&lambda;<sup>ae</sup><sub>im</sub>W<sup>&Lambda;</sup><sub>jebm</sub> + P(ij)P(ab)&lambda;<sup>a</sup><sub>i</sub>F<sup>&Lambda;</sup><sub>jb</sub> + P(ab)g<sub>ijae</sub>G<sub>be</sub> - + P(ij)g<sub>imab</sub>G<sub>mj</sub> 

**&Lambda;-CCSD Pseudoenergy**\
E<sub>&Lambda;</sub> = &lambda;<sup>a</sup><sub>i</sub>f<sub>ai</sub> + 0.25&lambda;<sup>ab</sup><sub>ij</sub>g<sub>abij</sub>

- - -
## CCSD -> CC2
CC2 is an approximate scheme to CCSD. The T<sub>1</sub> amplitudes are the same as CCSD, however the T<sub>2</sub> amplitudes are given by <&psi;<sup>ab</sup><sub>ij</sub>| **H<sub>N</sub>** e<sup>T<sub>1</sub></sup> + **F<sub>N</sub>** T<sub>2</sub>|0> = 0\
where **H<sub>N</sub>** and **F<sub>N</sub>** are the normal-ordered Hamiltonian and Fock matrices. **\[]** are the coupled-cluster diagram designations. Equations for T-amplitudes taken from reference above Table I(a) & (b) and (tilde) intermediates taken from Table III.

T<sub>1</sub> = f<sub>ai</sub> + *F*<sub>ae</sub>t<sup>e</sup><sub>i</sub> - *F*<sub>mi</sub>t<sup>a</sup><sub>m</sub> + *F*<sub>me</sub>t<sup>ae</sup><sub>im</sub> +
t<sup>e</sup><sub>m</sub>g<sub>amie</sub> - 0.5t<sup>ae</sup><sub>mn</sub>g<sub>mnie</sub> + 0.5t<sup>ef</sup><sub>im</sub>g<sub>amef</sub>
+ f<sub>ai</sub> **\[S<sub>1</sub>]**
+ *F*<sub>ae</sub>t<sup>e</sup><sub>i</sub>
+ + t<sup>e</sup><sub>i</sub>f<sub>ae</sub> **\[S<sub>3<sub>a</sub></sub>]**
+ + -0.5t<sup>e</sup><sub>i</sub>t<sup>a</sup><sub>m</sub>f<sub>me</sub> **\[S<sub>5<sub>a</sub></sub>]**
+ + t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>m</sub>g<sub>amef</sub> **\[S<sub>5<sub>b</sub></sub>]**
+ + -t<sup>e</sup><sub>i</sub> &tau;<sup>af</sup><sub>mn</sub>g<sub>mnef</sub>
+ + + -0.5t<sup>e</sup><sub>i</sub> t<sup>af</sup><sub>mn</sub>g<sub>mnef</sub> **\[S<sub>4<sub>a</sub></sub>]**
+ + + -0.5t<sup>e</sup><sub>i</sub> t<sup>a</sup><sub>m</sub>t<sup>f</sup><sub>n</sub>g<sub>mnef</sub> **\[S<sub>6</sub>]**
+ *F*<sub>mi</sub>t<sup>a</sup><sub>m</sub>
+ + -t<sup>a</sup><sub>m</sub>f<sub>mi</sub> **\[S<sub>3<sub>b</sub></sub>]**
+ + -0.5t<sup>a</sup><sub>m</sub>t<sup>e</sup><sub>i</sub>f<sub>me</sub> **\[S<sub>5<sub>a</sub></sub>]**
+ + -t<sup>a</sup><sub>m</sub>t<sup>e</sup><sub>n</sub>g<sub>mnie</sub> **\[S<sub>5<sub>c</sub></sub>]**
+ + -t<sup>a</sup><sub>m</sub> &tau;<sup>ef</sup><sub>in</sub>g<sub>mnef</sub>
+ + + -0.5t<sup>a</sup><sub>m</sub> t<sup>ef</sup><sub>in</sub>g<sub>mnef</sub> **\[S<sub>4<sub>b</sub></sub>]**
+ + + -0.5t<sup>a</sup><sub>m</sub> t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>n</sub>g<sub>mnef</sub> **\[S<sub>6</sub>]**
+ *F*<sub>me</sub>t<sup>ae</sup><sub>im</sub>
+ + t<sup>ae</sup><sub>im</sub>f<sub>me</sub> **\[S<sub>2<sub>a</sub></sub>]**
+ + t<sup>ae</sup><sub>im</sub>t<sup>f</sup><sub>n</sub>g<sub>mnef</sub> **\[S<sub>4<sub>c</sub></sub>]**
+ t<sup>e</sup><sub>m</sub>g<sub>amie</sub> **\[S<sub>3<sub>c</sub></sub>]**
+ -0.5t<sup>ae</sup><sub>mn</sub>g<sub>mnie</sub> **\[S<sub>2<sub>c</sub></sub>]**
+ 0.5t<sup>ef</sup><sub>im</sub>g<sub>amef</sub> **\[S<sub>2<sub>b</sub></sub>]**

T<sub>2</sub> = g<sub>abij</sub> + P(ij)(t<sup>ae</sup><sub>ij</sub>{*F*<sub>be</sub> - 0.5t<sup>b</sup><sub>n</sub>*F*<sub>me</sub>}) - P(ab)(t<sup>ab</sup><sub>im</sub>{*F*<sub>mj</sub> + 0.5t<sup>e</sup><sub>j</sub>*F*<sub>me</sub>}) + 0.5&tau;<sup>ab</sup><sub>mn</sub>*W*<sub>mnij</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub>*W*<sub>abef</sub> + P(ab)P(ij){t<sup>ae</sup><sub>im</sub>*W*<sub>mbej</sub> - t<sup>e</sup><sub>i</sub>t<sup>a</sup><sub>m</sub>g<sub>mbej</sub>} + P(ij)t<sup>e</sup><sub>i</sub>g<sub>abej</sub> - P(ab)t<sup>a</sup><sub>m</sub>g<sub>mbij</sub>

Since only T<sub>1</sub> operates on **H<sub>N</sub>** we can reduce the above equation to\
T<sub>2</sub> = g<sub>abij</sub> + 0.5&tau;<sup>ab</sup><sub>mn</sub>*W*<sub>mnij</sub> + 0.5&tau;<sup>ef</sup><sub>ij</sub>*W*<sub>abef</sub> - P(ab)P(ij)t<sup>e</sup><sub>i</sub>t<sup>a</sup><sub>m</sub>g<sub>mbej</sub> + P(ij)t<sup>e</sup><sub>i</sub>g<sub>abej</sub> - P(ab)t<sup>a</sup><sub>m</sub>g<sub>mbij</sub>
+ g<sub>abij</sub> **\[D<sub>1</sub>]**
+ 0.5&tau;<sup>ab</sup><sub>mn</sub>*W*<sub>mnij</sub> = t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>*W*<sub>mnij</sub> 
+ + t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>g<sub>mnij</sub> **\[D<sub>6<sub>b</sub></sub>]**
+ + t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>P(ij){t<sup>e</sup><sub>j</sub>g<sub>mnie</sub>} **\[D<sub>8<sub>b</sub></sub>]**
+ + 0.5t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>j</sub>g<sub>mnef</sub> **\[D<sub>9</sub>]**
+ 0.5&tau;<sup>ef</sup><sub>ij</sub>*W*<sub>abef</sub>
+ + t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>j</sub>g<sub>abef</sub> **\[D<sub>6<sub>a</sub></sub>]**
+ + -t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>j</sub>P(ab){t<sup>b</sup><sub>m</sub>g<sub>amef</sub>} **\[D<sub>8<sub>a</sub></sub>]**
+ + 0.5t<sup>e</sup><sub>i</sub>t<sup>f</sup><sub>j</sub>t<sup>a</sup><sub>m</sub>t<sup>b</sup><sub>n</sub>g<sub>mnef</sub> **\[D<sub>9</sub>]**
+ P(ab)P(ij)t<sup>e</sup><sub>i</sub>t<sup>a</sup><sub>m</sub>g<sub>mbej</sub> **\[D<sub>6<sub>c</sub></sub>]**
+ P(ij)t<sup>e</sup><sub>i</sub>g<sub>abej</sub> **\[D<sub>4<sub>a</sub></sub>]**
+ -P(ab)t<sup>a</sup><sub>m</sub>g<sub>mbij</sub> **\[D<sub>4<sub>b</sub></sub>]**

Need to add in the T<sub>2</sub> operating on **F<sub>N</sub>**
+ t<sup>ae</sup><sub>ij</sub>f<sub>be</sub> **\[D<sub>2<sub>a</sub></sub>]**
+ -t<sup>ab</sup><sub>im</sub>f<sub>mj</sub> **\[D<sub>2<sub>b</sub></sub>]**

- - -
## CCD and LCCD
The linear coupled-cluster doubles are a subset of the CCD equations. Hence there are no T<sub>1</sub> amplitudes and the T<sub>2</sub> amplitude are
+ g<sub>abij</sub> **\[D1]**
+ P(ab)*f*<sub>bc</sub>t<sup>ae</sup><sub>ij</sub> **\[D2a]**
+ -P(ij)*f*<sub>mj</sub>t<sup>ab</sup><sub>im</sub> **\[D2b]**
+ 0.5g<sub>abef</sub>t<sup>ef</sup><sub>ij</sub> **\[D2c]**
+ 0.5g<sub>mnij</sub>t<sup>ab</sup><sub>mn</sub> **\[D2d]**
+ P(ab)P(ij)g<sub>mbej</sub>t<sup>ae</sup><sub>im</sub> **\[D2e]**
+ 0.25g<sub>mnef</sub>t<sup>ef</sup><sub>ij</sub>t<sup>ab</sup><sub>mn</sub>  **\[D3a]**
+ P(ij)g<sub>mnef</sub>t<sup>ae</sup><sub>im</sub>t<sup>bf</sup><sub>jn</sub> **\[D3b]**
+ -0.5P(ij)g<sub>mnef</sub>t<sup>fe</sup><sub>im</sub>t<sup>ab</sup><sub>nj</sub> **\[D3c]**
+ -0.5P(ab)g<sub>mnef</sub>t<sup>ae</sup><sub>nm</sub>t<sup>fb</sup><sub>ij</sub> **\[D3d]**

For the Linear CCD we use **\[D1], \[D2a], \[D2b], \[D2c], \[D2d] and \[D2e]**
