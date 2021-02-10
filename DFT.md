# Density Functional Theory

Given that **harpy** contains much of the code needed for a DFT implementation I want to look into how it might be done.

For Kohm-Sham DFT we have \
E[&#x3C1;(**r**)] = T[&#x3C1;(**r**)] + V<sub>eN</sub>(&#x3C1;(**r**);R) +V<sub>NN</sub>(R) + V<sub>xc</sub>[&#x3C1;(**r**)], \
where &#x3C1;(**r**) is the density, T is kinetic energy, V<sub>eN</sub> is the electron-Nuclear attraction potential, V<sub>NN</sub> is the Nuclear-Nuclear repulsion potential and V<sub>xc</sub> is the exchange functional. 

We'll use LDA (Local Density Approximation) and the Vosko et al. xc-functional VWNIII (used for LDA in gaussian) \
```python
#exchange energy density
ex = -0.75 * pow(3.0 * rho /np.pi, 1/3)

#exchange potential - derivative wrt rho
Vx = - 4.0 * ex / 3.0

#radius of sphere with 1 electron at density rho
rs = 0.75 * pow(rho, -1/3) / np.pi

#fitting parameters
a = 0.0621814
b = 13.0720
c = 42.7198
xo = -0.409286

#parameters
x = np.sqrt(rs)
X = x*x + b*x + c
Q = np.sqrt(4.0*c - b*b)

#correlation energy
ec = 0.5*a*(np.log(x*x/X) + \
     2.0*b*math.atan(Q/(2.0*x+b)/Q - \
     b*x*(np.log((x-xo)*(x-xo)/X) + \
     2.0*(b+2.0*xo)*math.atan(Q/(2.0*x+b)/Q)/X
    
#correlation potential
vc = ec - a*(c*(x-xo) - b*x*x)/(6*((x-xo)*(x*x+b*x+c)))

```
Strategy \
1. Build basis
2. Overlap matrix (S)
3. Kinetic energy (K)
4. Nuclear attraction (J) 
5. Nuclear repulsion (ERI)
6. Core (J+K)

SG-1 radii \
```python
radii = {'H', 1.3; 'He', 'Li', 1.95; 'Be', 2.2 'B', 1.45; 'C', 1.2; \
         'N', 1.1; 'O', 1.1; 'F', 1.2; 'Na', 2.3; 'Mg', 2.3; \
         'Al', 2.1; 'Si', 1.3; 'P', 1.3; 'S', 1.1; 'Cl', 1.5;  }
 
lebedevPartitions = {['H', 'Li'] , [[6,6],[18,3],[26,1],[38,1],[74,1],[110,1],[146,6],[86,1],[50,1],[38,1],[18,1]] ; \
                     ['Be'], [[6,6],[18,2],[26,1],[38,2],[74,1],[86,1],[110,2],[146,5],[50,1],[38,1],[18,1],[6,2]] ; \
                     ['B'], [[6,4],[26,4],[38,3],[86,3],[146,6],[38,1],[6,2]] ; \
                     ['C'], [[6,6],[18,2],[26,1],[38,2],[50,2],[86,1],[110,1],[146,1],[170,2],[146,2],[86,1],[38,1],[18,1]] ; \
                     ['N'], [[6,6],[18,3],[26,1],[38,2],[74,2],[110,1],[170,2],[146,3],[86,1],[50,2]] ; \
                     ['O'], [[6,5],[18,1],[26,2],[38,1],[50,4],[86,1],[110,5],[86,1],[50,1],[38,1],[6,1]] ; \
                     ['F'], [[6,4],[38,2],[50,4],[74,2],[110,2],[146,2],[110,2],[86,3],[50,1],[6,1]] ; \
                     ['Na'], [[6,6],[18,2],[26,3],[38,1],[50,2],[110,8],[74,2],[6,2]] ; \
                     ['Mg'], [[6,5],[18,2],[26,2],[38,2],[50,2],[74,1],[110,2],[146,4],[110,1],[86,1],[38,2],[18,1],[6,1]] ; \
                     ['Al'], [[6,6],[18,2],[26,1],[38,2],[50,2],[74,1],[86,1],[146,2],[170,2],[110,2],[86,1],[74,1],[26,1],[18,1],[5,1]] ; \
                     ['Si', 'P'], [[6,5],[18,4],[38,4],[50,3],[74,1],[110,2],[146,1],[170,3],[86,1],[50,1],[6,1]] ; \
                     ['S'], [[6,4],[18,1],[26,8],[38,2],[50,1],[74,2],[110,1],[170,3],[146,1],[110,1],[50,1],[6,1]] ; \
                     ['Cl'], [[6,4],[18,7],[26,2],[38,2],[50,1],[74,1],[110,2],[170,3],[146,1],[110,1],[50,1],[6,1]] }
     
partitionNumbers = [23,23,23,23,23,23,23,23,26,26,26,26,26,26,26]
 
partitionCount = [1406,1406,1390,1426,1390,1414,1154,1494,1328,1492,1496,1496,1496,1456,1480]
 
#consistency checks
check = []
element = ['H','Li','Be','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl']
 
for s in element:
      
         
         
                      ```


