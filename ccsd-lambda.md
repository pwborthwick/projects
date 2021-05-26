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

