## Hugenholtz Diagrams

The idea is to analytically generate all the Hugenholtz diagrams of a given order. An approach is given [here](http://cds.cern.ch/record/609508/files/0303069.pdf). Also see Szabo and Ostlund (pg 256-262).

The perturbation series for the ground state energy can be given by diagrams in the Hugenholtz vein constructed according to the following rules 
For an *n*<sup>th</sup> order diagram, n dots (vertices) are drawn in a column. These vertices are connected by directed lines subject to the conditions
1. Each vertex has two lines pointing in and two pointing out.
2. Each diagram is connected, i.e.  one must be able to go from any one vertex to any other by following some number of lines.
3. No line connects a vertex with itself.
4. Each diagram is topologically distinct.

We will number vertices from 0 at the top, and down arrows are 'holes' (labelled a,b,c,...) and up arrows 'particles' (labelled r,s,t,...).

![](hugenholtz.png)

For order *n* there will be &half;*n*(*n*-1) pairs of nodes. So for order 2 there will be &half;2.1=1, and for 3 there will be &half;3.2=3 pairs as seen above.
```python
def nodalPairCount_hh(order):
    #return the number of pairs of nodes in Hugenholtz diagrams of 'order'
  
    return int(0.5 * order * (order- 1))
```   
For a diagram of order *n* there will be a total of 2*n* lines joining nodes in the diagram
```python
def nodalLineCount_hh(order):
    #return the total number of line connecting nodes in the diagram
    
    return 2 * order
```
First generate a list of nodal pairs, we will number the nodes 0, 1, 2, ...
```python
def nodalPairs_hh(order):
    #return number of nodal pairs in  diagram of 'order'
    
    pairs = []
    
    for p in range(order):
        for q in range(p+1, order):
            pairs.append([p, q])
            
    return pairs 
```
We need a list of lists of all combinations of connection between the pairs. Each element of the list will be the number of nodal pairs long. We generate it by a recursive subroutine which initially takes a list of zeros, nodal pairs long and a pair number 0 initially. So for order 3 it will take \[0,0,0] and 0. The routine then append to the initial list all possibilities on on the pair number ie \[0,0,0],\[1,0,0],\[2,0,0],\[3,0,0]. This list is then passed back to the routine with the pair number incremented, and so on until the pair number exceeds the number of nodal pairs. We can reduce the number of combinations by not allowing those where the sum of the connections exceeds the number of lines of the diagram order as any subsequent connections will still be above the limit. Once all pairs have been processed a routine is called to verify that the combination is valid according to the rules.
```python
connections = [[0] * nodalPairCount_hh(order)]
nodePair = 0
pairCount = nodalPairCount_hh(order)

def nodalPairConnectionsCombinations_hh(connections, nodePair, pairCount, order):
    #generate all combinations of connections between pairs
    
    #recursive return
    
    if nodePair == pairCount: 
        return validDiagrams_hh(connections,order)
    
    #define maximum number of connections
    limit = 3
    
    #adjust for special case order 2
    if order == 2: limit += 1
    
    #make copy of connections as we are modifying it and don't want to processed appended elements
    c = connections.copy()

    #loop over all elements in original connection list
    for connection in connections:    
    
        #loop over all possible connection types, 
        for i in range(1, limit+1):

            #make copy of current connection
            t = connection.copy()
            t[nodePair] = i

            #if sum of elements is less than or equal to allowed lines in diagram save
            if sum(t) <= 2 * order:
                c.append(t)
                
    #don't need original list now
    del connections
        
    #increment the nodePair
    nodePair += 1
      
    #recurse
    node = nodalPairConnectionsCombinations_hh(c, nodePair, pairCount, order)
        
    return node
```
Now need to check that connection combinations are valid. There are two checks 1. That the total number of connections equals the number of lines (2.order) and 2. That the diagram is connected ie you can get from node 0 to node order-1 via all other nodes. The first pair is always 0->1 so then first connection is from node 1. *tested* is an array initially set to False values that are set to true when a nodal pair has contributed to the route (the first element is True as we essentially start at node 1).
```python
def validDiagrams_hh(diagrams, order):
    #perform checks on diagrams for validity
    
    verified = []
    passed = False
    
    #get list of nodal pairs
    nodalPair = nodalPairs_hh(order)
    pairCount = nodalPairCount_hh(order)

    #number of lines
    for d in diagrams:
    
        if sum(d) == 2* order: 
        
            #check line count for diagram
            vertex = [0] * order
          
            #accummulate lines to each node
            for i in range(pairCount):
                vertex[nodalPair[i][0]] += d[i]
                vertex[nodalPair[i][1]] += d[i]
          
            #check correct number of lines at node
            if vertex == [4] * order:
          
                #number of lines at nodes verified check connected
                route = [0,1]
                node = 1
                tested = [False] * pairCount
                tested[0] = True
          
                #loop over nodes to be found
                while True:
          
                    #loop over nodal pairs
                    for i in range(pairCount):
            
                        pair = nodalPair[i]
              
                        #has pair been checked
                        if not tested[i]:
              
                            #is current node in nodal pair
                            if node == pair[0]: 
                                route.append(pair[1])
                                node = pair[1]
                            elif node == pair[1]: 
                                route.append(pair[0])
                                node = pair[0]
                            tested[i] = True
                                  
                    #all nodes tested leave while loop
                    if all(t == True for t in tested): break
                
                #have we got all nodes, make unique and sort
                route = list(set(route))
                route.sort()
                passed = True
                for i in range(order):
                    if i != route[i]: passed = False

                if passed: verified.append(d)
        
    return verified
```
We have our valid diagrams. Next determine the number of up arrow lines. The diagrams deemed valid are considered one at a time.
```python
def upArrowCombinations_hh(diagramCombinations,order):
    #get all combinations of up arrows

    for i in diagramCombinations:

        nodePair = -1
        pairCount = nodalPairCount_hh(order)

        arrows = [0] * pairCount
        upArrow_hh(arrows, nodePair, i, nodalPairs_hh(order), pairCount, order) 


    return
```
The diagram pairs are considered in turn and the combinations of possible up arrow configurations generated recursively. When all pairs have been processed the resulting up arrow configuration is tested for consistency.

```python
def upArrow_hh(up, nodePair, diagram, pairs, pairCount, order):
    #up arrow combinations for diagram

    if nodePair == pairCount-1: 

        passed = verifyArrow_hh(up, pairs, diagram, order)
        if passed: diagrams_hh.append([diagram,up.copy()])

        return 

    nodePair += 1

    #get limits of up connections
    lo = max(diagram[nodePair] - 2, 0)
    hi = min(diagram[nodePair],2)

    #generate combination within range
    for i in range(lo, hi+1):

        up[nodePair] = i
        upArrow_hh(up, nodePair, diagram, pairs, pairCount, order)

    nodePair -= 1

    return 
```
This routine tests for validity of an up arrow configuration. The number of up arrows at each node is determined and then checked for all nodes to have exactly 2 up arrow lines.
```python
def verifyArrow_hh(up, pairs, diagram, order):
    #check up arrow combination compatible with original diagram

    nodes = np.zeros(order)

    #loop over all pairs of nodes
    for n, pair in enumerate(pairs):
        i = pair[0]
        j = pair[1]

        #sum up arrows at each node
        if i < j:
            nodes[j] += up[n]
            nodes[i] += diagram[n] - up[n]
        else:
            nodes[i] += up[n]
            nodes[j] += diagram[n] - up[n]

    passed = True
    for i in range(order): 
        if nodes[i] != 2: passed = False 

    return passed
```
The whole procedure is run as
```python
import numpy as np

order = 4
diagrams_hh = []

pairCount = nodalPairCount_hh(order)
connections = [[0] * pairCount]
nodePair = 0

diagramCombinations = nodalPairConnectionsCombinations_hh(connections, nodePair, pairCount, order)
upArrowCombinations_hh(diagramCombinations, order)

d = []

for i, diagram in enumerate(diagrams_hh):
    if diagram[0]  != d:
        print('Base diagram ', diagram[0])
        d = diagram[0]  
    print('           arrow combination [',i,'] ', diagram[1])
```
for order 3...
```
Base diagram  [2, 2, 2]
           arrow combination [ 0 ]  [0, 2, 0]
           arrow combination [ 1 ]  [1, 1, 1]
           arrow combination [ 2 ]  [2, 0, 2]
```
for order 4...
```
Base diagram  [0, 1, 3, 3, 1, 0]
           arrow combination [ 1 ]  [0, 0, 2, 2, 0, 0]
           arrow combination [ 2 ]  [0, 1, 1, 1, 1, 0]
Base diagram  [0, 2, 2, 2, 2, 0]
           arrow combination [ 3 ]  [0, 0, 2, 2, 0, 0]
           arrow combination [ 4 ]  [0, 1, 1, 1, 1, 0]
           arrow combination [ 5 ]  [0, 2, 0, 0, 2, 0]
Base diagram  [0, 3, 1, 1, 3, 0]
           arrow combination [ 6 ]  [0, 1, 1, 1, 1, 0]
           arrow combination [ 7 ]  [0, 2, 0, 0, 2, 0]
Base diagram  [1, 0, 3, 3, 0, 1]
           arrow combination [ 8 ]  [0, 0, 2, 1, 0, 0]
           arrow combination [ 9 ]  [1, 0, 1, 2, 0, 1]
Base diagram  [2, 0, 2, 2, 0, 2]
           arrow combination [ 10 ]  [0, 0, 2, 0, 0, 0]
           arrow combination [ 11 ]  [1, 0, 1, 1, 0, 1]
           arrow combination [ 12 ]  [2, 0, 0, 2, 0, 2]
Base diagram  [3, 0, 1, 1, 0, 3]
           arrow combination [ 13 ]  [1, 0, 1, 0, 0, 1]
           arrow combination [ 14 ]  [2, 0, 0, 1, 0, 2]
Base diagram  [1, 3, 0, 0, 3, 1]
           arrow combination [ 15 ]  [0, 2, 0, 0, 1, 1]
           arrow combination [ 16 ]  [1, 1, 0, 0, 2, 0]
Base diagram  [2, 2, 0, 0, 2, 2]
           arrow combination [ 17 ]  [0, 2, 0, 0, 0, 2]
           arrow combination [ 18 ]  [1, 1, 0, 0, 1, 1]
           arrow combination [ 19 ]  [2, 0, 0, 0, 2, 0]
Base diagram  [3, 1, 0, 0, 1, 3]
           arrow combination [ 20 ]  [1, 1, 0, 0, 0, 2]
           arrow combination [ 21 ]  [2, 0, 0, 0, 1, 1]
Base diagram  [1, 1, 2, 2, 1, 1]
           arrow combination [ 22 ]  [0, 0, 2, 1, 0, 0]
           arrow combination [ 23 ]  [0, 1, 1, 0, 1, 0]
           arrow combination [ 24 ]  [0, 1, 1, 1, 0, 1]
           arrow combination [ 25 ]  [1, 0, 1, 1, 1, 0]
           arrow combination [ 26 ]  [1, 0, 1, 2, 0, 1]
           arrow combination [ 27 ]  [1, 1, 0, 1, 1, 1]
Base diagram  [1, 2, 1, 1, 2, 1]
           arrow combination [ 28 ]  [0, 1, 1, 0, 1, 0]
           arrow combination [ 29 ]  [0, 1, 1, 1, 0, 1]
           arrow combination [ 30 ]  [0, 2, 0, 0, 1, 1]
           arrow combination [ 31 ]  [1, 0, 1, 1, 1, 0]
           arrow combination [ 32 ]  [1, 1, 0, 0, 2, 0]
           arrow combination [ 33 ]  [1, 1, 0, 1, 1, 1]
Base diagram  [2, 1, 1, 1, 1, 2]
           arrow combination [ 34 ]  [0, 1, 1, 0, 0, 1]
           arrow combination [ 35 ]  [1, 0, 1, 0, 1, 0]
           arrow combination [ 36 ]  [1, 0, 1, 1, 0, 1]
           arrow combination [ 37 ]  [1, 1, 0, 0, 1, 1]
           arrow combination [ 38 ]  [1, 1, 0, 1, 0, 2]
           arrow combination [ 39 ]  [2, 0, 0, 1, 1, 1]
```
For example
```
Base diagram  \[3, 1, 0, 0, 1, 3]
           arrow combination \[ 20 ]  \[1, 1, 0, 0, 0, 2]           
           arrow combination \[ 21 ]  \[2, 0, 0, 0, 1, 1]
```           
corresponds to the two 4th order diagrams ![](310013.png)

To get the down arrows subtract the up arrows at each pair from the number of pairs at each node (take absolute value). So eg if \
\[3,1,0,0,1,3] is the base diagram and \[1,1,0,0,0,2] is the up arrow combination then \[2,0,0,0,1,1] is the down arrow combination. The code is
```python
def downArrow_hh(diagram, arrow, pairs):
    #compute the down arrows from an up arrow specification

    down = [0] * pairs
    for pair in range(pairs):
        down[pair] = abs(diagram[pair] - arrow[pair])
    
    return down
```
It is possible (but not easy) to draw the diagrams using pyplot and patches.Arc.

The first two rules for translating diagrams into algebraic expressions are:

1.  Each dot contributes an antisymmetric matrix element <in<sub>1</sub>in<sub>2</sub>||out<sub>1</sub>out<sub>2</sub>> to the numerator.
2.  Each pair of adjacent vertices contributes the denominator factor **&Sigma;** &epsilon;<sub>holes</sub> - **&Sigma;** &epsilon;<sub>particles</sub> 
    where the sums run over the labels of all hole and particle lines ceossing an imaginary line separating the pair of vertices.
    
In order to evaluate the above rules we need to know the nodes connecting a given node for both the in and out lines and the direction of the arrows on those lines. We will represent a join from node n to m as 'n:m' with n<m.
```python
def nodalFlows_hh(up, down, order):
    #collect the in and out lines at the node and their directions

    def sLine(node, p, d):
        #Construct line as 'node1node2' ordered

        join = [node, p]
        join.sort()
        return str(join[0]) + ':' + str(join[1]) + ':' + d

    l = []

    #loop over each node
    for node in range(order):

        _in   = []
        _out  = []

        #loop over pairs
        pairs = nodalPairs_hh(order)
        for i, pair in enumerate(pairs):

            #determine if in and out arrows are up or down
            if node == pair[0] and up[i] != 0:
                _in.append(sLine(node,pair[1], 'u'))
            if node == pair[1] and up[i] != 0:
                _out.append(sLine(node,pair[0], 'u'))
            if node == pair[0] and down[i] != 0:
                _out.append(sLine(node,pair[1], 'd'))
            if node == pair[1] and down[i] != 0:
                _in.append(sLine(node,pair[0], 'd'))


        l.append([_in, _out ])

    return l
```
This routine will return eg
```
flow at each node in - out
[[['0:1:u', '0:3:u'], ['0:1:d', '0:2:d']], [['0:1:d', '1:2:u'], ['0:1:u', '1:3:d']], 
[ ['0:2:d', '2:3:u'], ['1:2:u', '2:3:d']], [['1:3:d', '2:3:d'], ['0:3:u', '2:3:u']]]
```
From this the numerator can be read off directly as <'0:1:u''0:3:u'|'0:1:d''0:2:d'>... or substituting letters ('0:1:u'=r, '0:2:u' = s...'0:1:d'=a etc) then <rt||ab><au||re><bw|uf><ef||tw>. \
The denominators are found from
```python
def denominators_hh(connection, order):
    #Determine numerators of expression following rule 2

    def getNode(sNode,i):
        #decode line string 'n:m:d'
        if i == 0: return int(sNode[:sNode.find(':')])
        elif i == 1:
            idx = sNode.find(':')
            jdx = sNode.find(':', idx+1)
            return int(sNode[idx+1:jdx])
        elif i == 2:
            idx = sNode.find(':')
            jdx = sNode.find(':', idx+1)
            return sNode[jdx+1:]

    particles = []
    holes = []

    #loop over lines betwwen each pair of adjacent nodes
    for line in range(order-1):

        p = []
        h = []

        #enumerate of in and out connection pairs [[1:1,2:2],[3:3,4:4]]
        for i, pair in enumerate(connection):
            
            #enumerate each pair for in or out [1:1,2:2], [3:3,4:4]
            for j, arrow in enumerate(pair):

                #loop over joins in each pair 1:1, 2:2
                for join in range(2):

                    if getNode(arrow[join],0) <= line and getNode(arrow[join],1) > line:

                        if getNode(arrow[join],2) == 'u': p.append(arrow[join])
                        else: h.append(arrow[join])

        #use set to avoid double counting lines
        particles.append(list(set(p)))
        holes.append(list(set(h)))

    return particles, holes
```
which returns
```
Particles 
[['0:3:u', '0:1:u'], ['1:2:u', '0:3:u'], ['0:3:u', '2:3:u']]
holes
[['0:2:d', '0:1:d'], ['0:2:d', '1:3:d'], ['2:3:d', '1:3:d']]
```
or (**&epsilon;**<sub>b</sub>+**&epsilon;**<sub>a</sub>-**&epsilon;**<sub>t</sub>-**&epsilon;**<sub>r</sub>)(**&epsilon;**<sub>b</sub>+**&epsilon;**<sub>e</sub>-**&epsilon;**<sub>u</sub>-**&epsilon;**<sub>t</sub>)(**&epsilon;**<sub>f</sub>+**&epsilon;**<sub>e</sub>-**&epsilon;**<sub>t</sub>-**&epsilon;**<sub>w</sub>)

The next rule is \
The overall sign of the expression is (-1)<sup>h+l</sup> where h and l are the number of hole lines and closed loops, respectively. The closed loops are determined as follows: \
Write down the matrix elements eg 
<'0:1:u''0:3:u'|'0:1:d''0:2:d'><'0:1:d''1:2:u'|'0:1:u''1:3:d'><'0:2:d''2:3:u'|'1:2:u''2:3:d'><'1:3:d''2:3:d'|'0:3:u''2:3:u'> \
take element 1 **'0:1:u'** goto element 3 '0:1:d', goto next term element 2 '1:2:u' goto element 4 '1:3:d' continue with this pattern '0:2:d' -> '1:2:u' -> '2:3:d' -> '2:3:u' -> '0:3:u' -> '0:2:d' -> '0:1:d' -> **'0:1:u'**, when we get back to the starting element the number of elements in the sequence is the number of closed loops ie 11. The order is like <1.||..> -> <..||2.> -> <.3||..> -> <..||.4> -> <.5||..> -> <..||.6> *note when you cycle back to start use the other  number ie <.5||..> not <5.||..>. 
.


