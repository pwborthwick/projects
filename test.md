Need to create a program to test all aspects of harpy.</br>
[1] Need to ensure active molecule is what we need - for most part water sto-3g.\
Define geometry
```python
geometry = '''
              name=h2o
              matrix=c
              diis=off
              basis=sto-3g
              post={di}
              
              O1 8 0.000000000000 -0.143225816552 0.000000000000
              H1 1 1.638036840407 1.136548822547 -0.000000000000
              H2 1 -1.638036840407 1.136548822547 -0.000000000000

              end
           '''
 ```
 

[2] Replace harpy.hpf with geometry so copy harpy.hpf to harpy.tmp, create and write geometry to harpy.hpf. **remember to delete harpy.hpf and rename .tmp back to harpy.hpf**.</br>
```python
import os
os.rename('../source/harpy.hpf','harpy.hpf')
f = open('../source/harpy.hpf', 'w')
f.write(geometry)
f.close()
```
[3] Now we can create harpy.html with test results in it</br>
```python
import rhf
atoms,bases, data = rhf.mol(['geometry'])
e = rhf.scf(atoms, bases, data, ['postSCF'])
```
[4] Having created harpy.html we need to read it into a variable
```python
f = open('../source/harpy.html', 'r')
html = f.read()
```
Now we need to find keywords and associated values in the variable *html*. Examples are... </br>

     <table><tr><th>nuclear repulsion energy</th><th>0.71510434</th></tr></table>
     <table><tr><td>initial SCF electronic energy</td><td>-1.8380446</td></tr></table> 
     
These are relatively easy to parse</br>
But 

	  <table>
		  <caption>inertia tensor</caption>
		  <tr><th>x</th><th>y</th><th>z</th>
		  <tr><td>0.985423</td><td>0.0</td><td>0.0</td></tr>
		  <tr><td>0.0</td><td>0.985423</td><td>0.0</td></tr>
		  <tr><td>0.0</td><td>0.0</td><td>0.0</td></tr>
	  </table>
is more difficult.</br>
Ideas
```python
def getValues(html, key, element, n):
    #key is the keyword to search for, element is the delimiting HTML code and n is the number of values
    i = html.find(key)
    j = html[i:].find(element)
    k = html[i:].find(element[0]+'/'+element[1:])
    return html[j+4:k]
```    
[5] Things to check</br>
1. from geometry check principal momemts, rotational constants an nuclear repulsion energy.</br>
2. from e check electronic energy.</br>
3. check charges mu and ld, dipoles and quadrupoles.</br>
4.mp, om, ml, 




