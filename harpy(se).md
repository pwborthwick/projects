 ![image](https://user-images.githubusercontent.com/73105740/139957887-8532d90a-031c-4d12-940c-605cb21046b7.png)  <h2>Harpy (Semi-empirical)</h2>

### Why?
  When I was a post-graduate in the 1970's quantum chemistry was different than it is today. We had computers that needed the floor of a building to house them and they were slower than those available today. Because of this only limited calculations were possible by ab initio methods (Gaussian 70 was available) so in order to do calcuations, on the machines available, a plethora of semi-empirical methods arose. There was Extended Huckel Theory (EHT), Complete Neglect of Differential Overlap (CNDO), Intermediate Neglect of Differential Overlap (INDO), Neglect of Diatomic Differential Overlap (NDDO), Modified Intermediate Neglect of Differential Overlap (MINDO) and Perturbative Configuration Interaction using Localised Orbitals (PCILO) to name just a few. We used Pople and Beveridge's CNDO/2 program to do conformational studies using a CDC 7600 supercomputer. Of the semi-empirical methods CNDO and INDO were probably the most used and Gaussian still has an option to use the methods. Anyway I thought it would be fun to write a Python version of the original Pople and Beveridge program given in the appendix of 'Approximate Molecular Orbital Theory - JA Pople & DL Beveridge McGraw-Hill, 1970'. There are a few FORTRAN CNDO/INDO programs around eg [here](https://github.com/brhr-iwao/cindo_windows) but I haven't found a Python version.

### Installation
  Just copy all the files into a directory. You will need Python3, numpy and math but that's all. Once you've done that, from the installation directory type 'python validation.py' and you should see the test molecules each with a 'PASSED' by it's name.

  The molecule definitions are held in the file 'project.hpf'. Each molecule definition has the following form

     name=<string>
     [basis=slater]
   	 shell=closed|open
     [charge=<integer>]0
     [multiplicity=<integer>]1
     [tolerance=<float>]1e-8
     [cycles=<integer>]100
     [geometry=<string>]
     ]units=angstrom|bohr]bohr
     <blank line>
  Then number of atoms lines of atom identifier, atomic number, x-coordinate, y-coordinate, z-coordinate eg

     N1  7  1.651  0.000  0.000

  Followed by
  
     <blank line>
     end
     #------------------------------------------
     
   See file project.hpf for examples.

### Running
  To run the first molecule in project.hpf file use
     python 'harpy(se).py' -cndo
  or
     python 'harpy(se).py' -indo

  To run a specific molecule type eg
     python 'harpy(se).py' -LiH -cndo
  The molecule name is simply the 'name=' identifier in the project.hpf file. (if the identifier contains special characters enclose it in quotes eg -'oh-')

  To run from a different file use eg
     python 'harpy(se).py' -water.hpf -indo
  The file must have a .hpf extension to be recognized.

### Files
  The files included are/
  **harpy(se).py** - handles the command line arguements and passes them to cindo.py, also does the timing of the run and prints it and returned energy to console.

  **cindo.py** - project.hpf commands parsing, overlap and coulomb integral evaluation, calls closed.py or open.py routines.

  **integral.py** - reduced overlap evaluation and auxillary functions.

  **atom.py** - atom class and various helper functions, harmonic transform matrix.

  **data.py** - the semi-emirical values are held here, output data formatting and printing.

  **closed.py** - closed shell Huckel Hamiltonian, closed shell scf iterations.

  **open.py** - open shell Huckel Hamiltonian, open shell scf iterations.

  **validation.py** - runs a suite of test molecules provided in project.hpf. Compares total energy and orbital eigenvectors against known correct values. This is a test of the self-consistency of the program not of the program against known CNDO/2 or INDO results. Although see below.
  
  **project.hpf** - this holds the molecular geometries and parameters for each molecule calculation..

### Output 
  The output (to console) contains
+ Echo of input file parameters
+ List of atoms, their coordinates (bohr), orbital number range assigned to atom, Slater parameters
+ Optional requested geometry<sup>\*</sup>
+ List of orbitals, their atom center, quantum numbers, designation (eg s, dxy), 'occupied' or ''
+ Overlap matrix
+ Coulomb matrix
+ Core Hamiltonian matrix
+ SCF interations, iteration number, electronic energy (Hartrees), delta energy, trend<sup>\*\*</sup>
+ Final electronic, nuclear, total and binding energies
+ Eigen values and eigenvectors
+ Density matrix
+ Atomic populations and charges
+ Dipole moments from density, sp, pd and total components
+ Inferred bonds and Wiberg-Mayer bond-order
open shell
+ &alpha; and &beta; electron densities and charges
+ spin density matrix
+ &alpha;  and &beta; spin densities and hyperfine coupling constants

<sup>\*</sup> The geometry specifier in the project file is of the form 'geometry=\<specifier>...\<specifier>', where \<specifier> can be either\
	{r:i,j} or {a:i,j,k} or {t:i,j,k,l}, where i,j,k,l are integers denoting atoms (first set of coordinates is atom 1). 'r' is bond, 'a' is angle and 't' is torsion angle.

<sup>\*\*</sup> Trend is a string of the sign of the difference of the last two scf electronic energies. If all is well it should be a sequence like '-','--','---' etc. 

### Testing
  I have checked the results against output files contained in this [download](http://www.jh-inst.cas.cz/~liska/Cnindo.htm). (The input and output files are a mess often the contents aren't what the file title advertises!) There is also an on-line CNDO/2 application [here](https://www.colby.edu/chemistry/PChem/scripts/cndo.html). Harpy(se) has been tested against these sources. The verification.py program contains reference values from these programs. If you are are looking at these programs be aware that with the 'jh-inst' program if the iteration limit is exceeded it just uses the unconverged value to print post-scf results without any warning - there's just no 'energy satisfied' message. The colby.edu program has a very low convergence tolerance (1e-4?) so results may not be too accurate.
