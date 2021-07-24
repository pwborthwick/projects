## VQE Using Sympy
  Josh Goings in his blog wrote about the [Variational Quantum Eigensolver](https://joshuagoings.com/blog/) based on [this paper](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.6.031007). I wanted to see how far I could get replicating Josh's article using sympy.quantum.
+ **qubits**
 Single qubit basis states are written |0> and |1> and have a vector representation of 
 
 ![image](https://user-images.githubusercontent.com/73105740/126635519-d249b848-b4b7-4e3f-bc39-2df70f6061bf.png).
 
 In Sympy we can write
 ```python
>>>_0 = Qubit('0')
>>>print(_0,qubit_to_matrix(_0))
|0> Matrix([[1], [0]])

>>>_1 = Qubit('1')
>>>print(_1,qubit_to_matrix(_1))
|1> Matrix([[0], [1]])
```
Josh shows that the outer product of basis qubits are

![image](https://user-images.githubusercontent.com/73105740/126637289-d40f946a-9d37-40bb-b25b-f916fa7fc7f1.png)

 In Sympy we can verify this, first we have to get bras from our kets
 ```python
>>>_0t = Dagger(_0)
>>>print(_0t,qubit_to_matrix(_0t))
<0| Matrix([[1, 0]])
```
Now we can compute the outer product
```python
>>>_0t = Dagger(_0)
>>>_00t = _0 * _0t
>>>print(_00t,qubit_to_matrix(_00t))

>>>_1t = Dagger(_1)
>>>print(_0t,qubit_to_matrix(_0t))
>>>_11t = _1*_1t
>>>print(_11t,qubit_to_matrix(_11t))
|0><0| Matrix([[1, 0], [0, 0]])
|1><1| Matrix([[0, 0], [0, 1]])
```
+ **States**
 Josh writes a single qubit state composed of 2 basis qubits
 
 ![image](https://user-images.githubusercontent.com/73105740/126638881-2d6cc2e0-de22-4e79-a192-7fe3501ca3a0.png)
 
 We write this as 
 ```python
 >>>psi = Ket(_0, _1)
 >>>print(psi)
 ||0>|1>>
 ```
+ **Gates**
 Gates are matrices that operate on single qubit states. Josh shows
 
 ![image](https://user-images.githubusercontent.com/73105740/126640071-1753aef0-2ffe-4aeb-8a46-d940d83dbbd7.png)

Below the (1) arguments are the number of qubits the gate operates on, here single qubits
```python
>>>print('Hadamard Gate ', H.get_target_matrix(H))
>>>print('Identity Gate ', IdentityGate.get_target_matrix(1))
>>>print('Pauli X ',XGate.get_target_matrix(1))
>>>print('Pauli Y', YGate.get_target_matrix(1))
>>>print('Pauli Z ',ZGate.get_target_matrix(1))
>>>print('Phase Gate ',S.get_target_matrix(1))

Hadamard Gate  Matrix([[1/sqrt(2), 1/sqrt(2)], [1/sqrt(2), -sqrt(2)/2]])
Identity Gate  Matrix([[1, 0], [0, 1]])
Pauly X  Matrix([[0, 1], [1, 0]])
Pauli Y Matrix([[0, -I], [I, 0]])
Pauli Z  Matrix([[1, 0], [0, -1]])
Phase Gate  Matrix([[1, 0], [0, I]])
```
Josh also shows the rotation matrices

![image](https://user-images.githubusercontent.com/73105740/126642576-10db9161-b3af-4321-9506-435be2f45fb2.png)

There doesn't seem to be a built in version of these but we can construct them
```python
>>>theta = Symbol('theta')

>>>Rx = Matrix([[cos(theta/2), -I*sin(theta/2)],[-I*sin(theta/2),cos(theta/2)]])
>>>Ry = Matrix([[cos(theta/2), -sin(theta/2)],[sin(theta/2),cos(theta/2)]])
>>>Rz = Matrix([[exp(-I*theta/2), 0],[0,exp(I*theta/2)]])

>>>print(Rx)
>>>print(Ry)
>>>print(Rz)

Matrix([[cos(theta/2), -I*sin(theta/2)], [-I*sin(theta/2), cos(theta/2)]])
Matrix([[cos(theta/2), -sin(theta/2)], [sin(theta/2), cos(theta/2)]])
Matrix([[exp(-I*theta/2), 0], [0, exp(I*theta/2)]])
```
+ **2-qubit states**
 Josh goes on to give
 
 ![image](https://user-images.githubusercontent.com/73105740/126643620-db9d3ccd-a1f4-4e36-8410-c2987621383a.png)
 
 We implement those as 
 ```python
>>>print(Qubit('00'), qubit_to_matrix(Qubit('00')))
>>>print(Qubit('01'), qubit_to_matrix(Qubit('01')))
>>>print(Qubit('10'), qubit_to_matrix(Qubit('10')))
>>>print(Qubit('11'), qubit_to_matrix(Qubit('11')))

|00> Matrix([[1], [0], [0], [0]])
|01> Matrix([[0], [1], [0], [0]])
|10> Matrix([[0], [0], [1], [0]])
|11> Matrix([[0], [0], [0], [1]])
```
We can do the kronecker product way
```python
>>>print(qubit_to_matrix(TensorProduct(Qubit('0'),Qubit('0'))))
Matrix([[1], [0], [0], [0]])
```
+ **More Gates**
 Josh now introduces the CNOT and SWAP gates. CNOT<sub>01</sub>, the first subscript is the 'control' qubit and the second is the 'target' qubit. Depending on the value of the control qubit the target qubit is changed. Consider the action of CNOT<sub>01</sub> on each of <00|, <01|, <10| and <11|
 ```python
>>>print('00',qapply(CNotGate(1,0)*Qubit('00')))
>>>print('01',qapply(CNotGate(1,0)*Qubit('01')))
>>>print('10',qapply(CNotGate(1,0)*Qubit('10')))
>>>print('11',qapply(CNotGate(1,0)*Qubit('11')))

00 |00>
01 |01>
10 |11>
11 |10>
```
So if the first (control) qubit is 0 nothing happens to the second (target) qubit, but if the first qubit is 1 the second qubit is changes. There is also CNOT<sub>10</sub>. Swap is defined
```python
>>>print('Swap ',SwapGate.get_target_matrix(SwapGate))

Swap  Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
```
This takes |00> -> |00>, |01> -> |10>, |10> -> |10> and |11> -> |11>\
There is an example given

![image](https://user-images.githubusercontent.com/73105740/126678933-fec5109c-c755-4f26-8766-2acf8d0900b3.png)

We can check this
```python
>>>print(TensorProduct(XGate.get_target_matrix(2),YGate.get_target_matrix(2)))

Matrix([[0, 0, 0, -I], [0, 0, I, 0], [0, -I, 0, 0], [I, 0, 0, 0]])
```
+ **Classical Result**
 The Hamiltonian given is (see Josh's blog for details)
 
 ![image](https://user-images.githubusercontent.com/73105740/126679742-980c3e8c-db48-417b-95eb-050ad1774f41.png)

This is implemented as
```python
>>>g0, g1, g2, g3, g4, g5 = [-0.4804, +0.3435, -0.4347, +0.5716, +0.0910, +0.0910]

>>>ti = IdentityGate.get_target_matrix(2)
>>>t0 = TensorProduct(ti,ti)
>>>tz = ZGate.get_target_matrix(2)
>>>t1 = TensorProduct(ti, tz)
>>>t2 = TensorProduct(tz, ti)
>>>t3 = TensorProduct(tz, tz)
>>>ty = YGate.get_target_matrix(2)
>>>t4 = TensorProduct(ty, ty)
>>>tx = XGate.get_target_matrix(2)
>>>t5 = TensorProduct(tx, tx)
>>>H = g0*t0 + g1*t1 + g2*t2 + g3*t3 + g4*t4 + g5*t5
>>>print(H)

>>>nuclear_repulsion = 0.7055696146

>>>electronic_energy = H.eigenvals()
>>>vals = list(electronic_energy.keys())
>>>vals.sort()

>>>print("Classical diagonalization: {:+2.8} Eh".format(vals[0] + nuclear_repulsion))
>>>print("Exact (from G16):          {:+2.8} Eh".format(-1.1457416808))

[1.11022302462516e-16, 0,                  0,                 0                ], 
[0,                   -1.83020000000000,   0.182000000000000, 0                ], 
[0,                    0.182000000000000, -0.273800000000000, 0                ], 
[0,                    0,                  0,                 0.182400000000000]])
Classical diagonalization: -1.1456295 Eh
Exact (from G16):          -1.1457417 Eh
```
This is the same as Josh's value (not surprisingly as we've used the same values) but this so far is not using qubits on a quantum circuit.

+ **Initial State**
 The initial state in this problem with just 2 electrons is |01>. We can generate this directly as Qubit('01') but Josh derives this from <00| by 1) Defining a zero 4-vector, 2) setting position 1 to 1. 3) acting on that vector with tensor product of identity and S<sub>x</sub>.
 ```python
>>>print('Initial State ', Qubit('01'))
>>>print(TensorProduct(ti,tx)*Matrix([1,0,0,0]))
Initial State  |01>
Matrix([[0], [1], [0], [0]])
```
+ **Ansatz**
 This problem is tackled using the Unitary Coupled Cluster (UCC) ansatz, viz. (in this example U(&theta;) = exp(-i&theta;X<sub>0</sub>Y<sub>1</sub>) where X<sub>0</sub> is X-gate acting on the 0 (top) qubit. The parameterized H<sub>2</sub> wave function is then |&psi;(&theta;)> = exp(-i&theta;X<sub>0</sub>Y<sub>1</sub>) |01>\
 ```python
 ansatz = lambda theta: exp(-I* theta*TensorProduct(ty,tx))
```
 The expected value of the Hamiltonian is then <&psi;|H|&psi;> or **&psi;**<sup>T</sup>**H&psi;**
```python
def expected(theta,ansatz,Hmol,psi0):
    circuit = ansatz(theta[0])
    psi = circuit*psi0
    return (psi.transpose().conjugate() * (H*psi))[0,0]
```
The whole of Josh's 'lazy VQE' in Sympy
```python
>>>def expected(theta,ansatz,Hmol,psi0):
>>>    circuit = ansatz(theta[0])
>>>    psi = circuit*psi0
>>>    return (psi.transpose().conjugate() * (H*psi))[0,0]

>>>ansatz = lambda theta: exp(-I* theta*TensorProduct(ty,tx))
>>>from scipy.optimize import minimize

>>># initial guess for theta
>>>theta  = [0.0]
>>>result = minimize(expected,theta,args=(ansatz,H,psi0))
>>>theta  = result.x[0]
>>>val    = result.fun

>>>print("Lazy VQE: ")
>>>print("  [+] theta:  {:+2.8} deg".format(theta))
>>>print("  [+] energy: {:+2.8} Eh".format(val + nuclear_repulsion))
Lazy VQE: 
  [+] theta:  -0.11487186 deg
  [+] energy: -1.1456295 Eh
```
+ **Quantum Ansatz**
 


 
