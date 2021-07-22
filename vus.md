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

```python
>>>print('Hadamard Gate ', H.get_target_matrix(H))
>>>print('Identity Gate ', IdentityGate.get_target_matrix(1))
>>>print('Pauly X ',XGate.get_target_matrix(XGate))
>>>print('Pauli Y', YGate.get_target_matrix(YGate))
>>>print('Pauli Z ',ZGate.get_target_matrix(ZGate))
>>>print('Phase Gate ',S.get_target_matrix(S))

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





 
