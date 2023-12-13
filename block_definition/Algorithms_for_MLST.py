
# coding: utf-8

# In[58]:

import queue
import matplotlib.pyplot as plt
import pdb


# In[59]:

#Creates a small example adjacency list.
def exampleGraph1():
    fullgraph=dict()
    fullgraph[1]=[-1,2,-1,3,-1,-1]
    fullgraph[2]=[1,-1,-1,4,-1,-1]
    fullgraph[3]=[-1,4,1,7,-1,5]
    fullgraph[4]=[3,-1,2,8,-1,6]
    fullgraph[5]=[-1,6,-1,9,3,11]
    fullgraph[6]=[5,16,-1,10,4,12]
    fullgraph[7]=[-1,8,3,-1,-1,9]
    fullgraph[8]=[7,-1,4,-1,-1,10]
    fullgraph[9]=[-1,10,5,-1,7,13]
    fullgraph[10]=[9,15,6,-1,8,14]
    fullgraph[11]=[-1,12,-1,13,5,-1]
    fullgraph[12]=[11,-1,-1,14,6,-1]
    fullgraph[13]=[-1,14,11,-1,9,-1]
    fullgraph[14]=[13,-1,12,-1,10,-1]
    fullgraph[15]=[10,17,16,23,-1,21]
    fullgraph[16]=[6,18,-1,15,-1,19]
    fullgraph[17]=[15,-1,18,24,-1,22]
    fullgraph[18]=[16,-1,-1,17,-1,20]
    fullgraph[19]=[-1,20,-1,21,16,-1]
    fullgraph[20]=[19,-1,-1,22,18,-1]
    fullgraph[21]=[-1,22,19,-1,15,-1]
    fullgraph[22]=[21,-1,20,25,17,-1]
    fullgraph[23]=[-1,24,15,-1,-1,-1]
    fullgraph[24]=[23,-1,17,-1,-1,25]
    fullgraph[25]=[-1,-1,22,-1,24,-1]
    return fullgraph


# In[60]:

def exampleGraphSymmetry():
    d=dict()
    d[0]=[37,19,28,1,-1,10]
    d[1]=[-1,-1,0,2,-1,-1]
    d[2]=[-1,-1,1,4,-1,-1]
    d[3]=[-1,4,-1,-1,-1,-1]
    d[4]=[3,6,2,-1,-1,5]
    d[5]=[-1,-1,-1,-1,4,-1]
    d[6]=[4,8,-1,7,-1,-1]
    d[7]=[-1,-1,6,-1,-1,-1]
    d[8]=[6,-1,-1,-1,9,-1]
    d[9]=[-1,-1,-1,-1,-1,8]
    d[10]=[-1,-1,-1,-1,0,11]
    d[11]=[-1,-1,-1,-1,10,12]
    d[12]=[13,15,14,-1,11,-1]
    d[13]=[-1,12,-1,-1,-1,-1]
    d[14]=[-1,-1,-1,12,-1,-1]
    d[15]=[12,17,-1,-1,-1,16]
    d[16]=[-1,-1,-1,-1,15,-1]
    d[17]=[15,-1,-1,18,-1,-1]
    d[18]=[-1,-1,17,-1,-1,-1]
    d[19]=[0,20,-1,-1,-1,-1]
    d[20]=[19,24,-1,-1,-1,-1]
    d[21]=[-1,-1,22,-1,-1,-1]
    d[22]=[-1,-1,-1,21,-1,23]
    d[23]=[-1,27,-1,-1,22,24]
    d[24]=[20,-1,26,-1,23,25]
    d[25]=[-1,-1,-1,-1,24,-1]
    d[26]=[-1,-1,-1,24,-1,-1]
    d[27]=[23,-1,-1,-1,-1,-1]
    d[28]=[-1,-1,29,0,-1,-1]
    d[29]=[-1,-1,35,28,-1,-1]
    d[30]=[-1,31,-1,-1,-1,-1]
    d[31]=[30,-1,-1,-1,32,-1]
    d[32]=[-1,-1,34,-1,35,31]
    d[33]=[35,-1,-1,-1,-1,-1]
    d[34]=[-1,-1,-1,32,-1,-1]
    d[35]=[-1,33,-1,29,36,32]
    d[36]=[-1,-1,-1,-1,-1,35]
    d[37]=[-1,0,-1,38,-1,-1]
    d[38]=[-1,-1,37,-1,-1,-1]
    return d


# In[61]:

def bfs(fullgraph,startnode):
    '''
    Runs BFS on a graph from a given start node.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/19/2017
    
    Keyword arguments:
    fullgraph -- a graph in adjacency list representation stored as a dictionary
    startnode -- the dictionary key of the node that we will start at

    Executes BFS algorithm, returning the max distance from the startnode to any other node, and the associated tree.
    '''
    dist=dict()
    #infinite distance (size of the dictionary works because the 
    #distance between two nodes is strictly less than the number of nodes)
    inf=len(fullgraph)
    for key in fullgraph:
        dist[key]=inf
    dist[startnode]=0
    queue1=queue.Queue(0)
    queue1.put(startnode)
    current=0
    maxdist=0
    maxkey=[]
    while not queue1.empty():
        current=queue1.get()
        for adjacent in fullgraph[current]:
            if adjacent!=-1 and dist[adjacent]==inf:
                queue1.put(adjacent)
                dist[adjacent]=dist[current]+1
    for key in dist:
        if dist[key]>maxdist:
            maxdist=dist[key]
            maxkey=[key]
        elif dist[key]==maxdist:
            maxkey.append(key)
    return {'maxdist': maxdist, 'maxkey': maxkey, 'distances' : dist}


# In[62]:

class Pod:
    '''
    Helper class for DLList class. Creates pods that hold data and can be connected to each other.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/20/2017
    
    Class that creates Pod objects. Pod objects remember a value, and 2 adjacent pods (previous,
    and subsequent pods). Used in the doubly linked list implemented in DLList class.
    '''
    def __init__(self, previous, subsequent, value):
        self.previous=previous
        self.subsequent=subsequent
        self.value=value
    def setPrevious(self, previous):
        self.previous=previous
    def getPrevious(self):
        return self.previous
    def setSubsequent(self, subsequent):
        self.subsequent=subsequent
    def getSubsequent(self):
        return self.subsequent
    def setValue(self, value):
        self.setValue()
    def getValue(self):
        return self.value


# In[63]:

#Class that creates a doubly linked list data structure in support of the mslt (Maximum Leaf Spanning Tree)
#function. Essentially, provides a framework for the instances of the Pod class to work with.
#Remembers the first and last pods in the list.
class DLList:
    '''
    Class that creates and manages a doubly linked list.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/20/2017
    
    Class that creates a doubly linked list. The list is formed of Pod objects that are connected and
    assigned values. A DLList instance remembers the first and las pod in its list.
    '''
    
    def __init__(self):
        '''
        Constructor. List starts out empty, so there is no first or last pod.

        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

        Date: 7/20/2017
        '''
        self.first=None
        self.last=None
    
    def __str__(self):
        '''
        Returns a string representation of the list.

        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

        Date: 7/20/2017
        '''
        if self.empty():
            return "<>"
        output="<"
        current=None
        nextpod=self.first
        while nextpod!=None:
            current=nextpod
            nextpod=current.getSubsequent()
            output=output+str(current.getValue())+", "
        output=output[:-2]+">"
        return output
        
    def empty(self):
        '''
        Checks to see if the list is empty.

        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

        Date: 7/20/2017
        '''
        if self.first==None:
            return True
        return False
    
    def add(self, value):
        '''
        Adds a value to the end of the list.

        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

        Date: 7/20/2017
        
        Keyword arguments:
        value -- The value to be added to the list.
        
        Adds a value to the end of the list by putting it in a new pod and attaching this pod
        to the previous last pod.
        '''
        temp=Pod(self.last,None,value)
        if self.first == None:
            self.first=temp
            self.last=temp
        else:
            self.last.setSubsequent(temp)
            self.last=temp
        return temp
    
    def pop(self):
        '''
        Returns the value at the front of the list and removes the pod containing it from the list

        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

        Date: 7/20/2017
        '''
        temp=self.first
        if temp is not self.last:
            self.first=temp.getSubsequent()
            self.first.setPrevious(None)
        else:
            self.first=None
            self.last=None
        return temp.getValue()
    
    def deletePod(self, current):
        '''
        Deletes a specified pod from the list.

        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

        Date: 7/20/2017
        
        Keyword arguments:
        current -- The pod to be deleted.
        
        Deletes a specified pod from the list. "current" points to the pod that needs to be deleted,
        so deletion thus consists of connecting the precedecing and following pods.
        '''
        if current is self.first and current is self.last:
            self.first=None
            self.first=None
        elif current is self.first:
            self.first=current.getSubsequent()
            self.first.setPrevious(None)
        elif current is self.last:
            self.last=current.getPrevious()
            self.last.setSubsequent(None)
        else:
            current.getPrevious().setSubsequent(current.getSubsequent())
            current.getSubsequent().setPrevious(current.getPrevious())


# In[64]:

def expand(vertex, forest, inforest, fullgraph, outsidedegrees, pointers, whichrule, rule1, rule2, rule3, rule4, numFaces):
    '''
    Helper method for mlst, that expands and updates rule lists.

    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

    Date: 7/20/2017
        
    Keyword arguments:
    vertex -- The vertex we wish to expand
    forest -- The forest we are building
    inforest -- Dictionary that keeps track of whether each vertex is in our forest
    fullgraph -- The original graph
    outsidedegrees -- The number of neighbors that a vertex has that are not in our forest
    pointers -- Dictionary that keeps track of which pod is associated with a vertex
    whichrule -- Dictionary that keeps track of which rule is applicable to a vertex
    rule1 -- DLList of vertices that qualify for rule 1
    rule2 -- DLList of vertices that qualify for rule 2
    rule3 -- DLList of vertices that qualify for rule 3
    rule4 -- DLList of vertices that qualify for rule 4
        
    Helper method for mlst. Implements the expand method described in Solis-Oba et al, as well
    as implements the process of updating the rule lists described in same.
    '''
    #we will keep track of the vertices that had their outside degrees (number of immediate
    #neighbors not in the forest) change or were added to the forest, because we will need to
    #update our rule lists to reflect these changes
    #print(vertex)
    #print(outsidedegrees)
    #if vertex==17:
    #    print(rule1)
    #    print(whichrule[18])
    degreechanged=set()
    addedtoforest=[]
    pointers[vertex]=None
    connected=False
    #if the vertex is not in our forest, we add it to our forest and decrement the outside degrees
    #of its neighbors.
    if not inforest[vertex]:
        inforest[vertex]=True
        addedtoforest.append(vertex)
        for neighbor2 in fullgraph[vertex]:
            if neighbor2!=-1:
                outsidedegrees[neighbor2]=outsidedegrees[neighbor2]-1
                degreechanged.add(neighbor2)
    #we add the neighbors to the forest if they are not in the forest and adjust outside degrees.
    for i in range(numFaces):
        neighbor=fullgraph[vertex][i]
        if neighbor!=-1:
            if not inforest[neighbor]:
                inforest[neighbor]=True
                forest[vertex][i]=neighbor
                if numFaces==6:
                    forest[neighbor][(i%2+1)%2+i-(i%2)]=vertex
                else:
                    forest[neighbor][fullgraph[neighbor].index(vertex)]=vertex
                addedtoforest.append(neighbor)
                for neighbor2 in fullgraph[neighbor]:
                    if neighbor2!=-1:
                        outsidedegrees[neighbor2]=outsidedegrees[neighbor2]-1
                        degreechanged.add(neighbor2)
            #if we are dealing with a rule 2 or a rule 3 expansion, we need to connect our vertex to one already
            #in the forest. So we take the first neighbor we find that is in the forest, and connect it
            elif (whichrule[vertex] == 2 or whichrule[vertex] == 3) and not connected:
                connected=True
                forest[vertex][i]=neighbor
                if numFaces==6:
                    forest[neighbor][(i%2+1)%2+i-(i%2)]=vertex
                else:
                    forest[neighbor][fullgraph[neighbor].index(vertex)]=vertex
    #We update our rule lists. If it is now in the forest, we must remove it from 2,3 or 4.
    #If it satisfies the requirements for 1, then we add it to 1
    whichrule[vertex]=0
    for key in addedtoforest:
        if key != vertex:
            if whichrule[key] == 2:
                rule2.deletePod(pointers[key])
            elif whichrule[key] == 3:
                rule3.deletePod(pointers[key])
            elif whichrule[key] == 4:
                rule4.deletePod(pointers[key])
        if outsidedegrees[key] >= 2:
            whichrule[key]=1
            pointers[key]=rule1.add(key)
        else:
            whichrule[key]=0
            pointers[key]=None
    #We continue to update our rule lists. If its degree changed then we make sure the vertex is in the correct list.
    for key in degreechanged:
        if not inforest[key]:
            if whichrule[key] == 2:
                rule2.deletePod(pointers[key])
            elif whichrule[key] == 3:
                rule3.deletePod(pointers[key])
            elif whichrule[key] == 4:
                rule4.deletePod(pointers[key])
            if outsidedegrees[key] >= 3:
                whichrule[key]=2
                pointers[key]=rule2.add(key)
            elif outsidedegrees[key] == 2:
                whichrule[key]=3
                pointers[key]=rule3.add(key)
            else:
                whichrule[key]=0
                pointers[key]=None
        elif key!=vertex:
            #updating rule1 list is not described in the paper, however updating it seems necessary.
            if whichrule[key] == 1 and outsidedegrees[key]<2:
                #print(key)
                #print(vertex)
                #print(rule1)
                rule1.deletePod(pointers[key])
                whichrule[key]=0
                pointers[key]=None
    #print(forest)
    #print(whichrule)
    #print("rule1:",rule1)
    #print("rule4:",rule4)


# In[65]:

def mlst(fullgraph, numFaces):
    '''
    A linear 2-approximation algorithm for approximating a MSLT.

    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

    Date: 7/20/2017
        
    Keyword arguments:
    fullgraph - A dictionary that is an adjacency list for a graph of possible connections
        
    A linear 2-approximation algorithm for approximating a MSLT (Maximum Leaf Spanning Tree) as
    described in Solis-Oba et al.
    '''
    #the forest we will be building (in adjacency list format)
    forest=dict()
    #lists of vertices that qualify for one of the 4 rules described in the paper. The lists are the doubly
    #linked lists described above.
    rule1=DLList()
    rule2=DLList()
    rule3=DLList()
    rule4=DLList()
    #dictionary that stores the number representing which rule can be applied to each vertex (0 for no rule)
    whichrule=dict()
    #dictionary that contains the Pods (described above) for each vertex's place in one of the rule lists.
    #if the vertex does not qualify for any rules, it has None instead of a Pod.
    pointers=dict()
    #the number of immediate neighbors a vertex has that are not in our forest
    outsidedegrees=dict()
    #dictionary containing booleans of whether or not a vertex is in our forest.
    inforest=dict()
    #sets up all our dictionaries & lists, and figures out which vertices qualify to start in the
    #rule 4 lists (it is impossible for vertices to qualify for other lists because our forest is empty)
    for key in fullgraph:
        forest[key]=[-1]*numFaces
        inforest[key]=False
        outsidedegrees[key]=0
        for i in range(numFaces):
            if fullgraph[key][i]!=-1:
                outsidedegrees[key]+=1
        if outsidedegrees[key] >= 3:
            pointers[key]=rule4.add(key)
            whichrule[key]=4
        else:
            pointers[key]=None
            whichrule[key]=0
    if rule4.empty():
        return dict()
    #While we can expand our forest, we do so according to the rules described in the paper.
    while not (rule1.empty() and rule2.empty() and rule3.empty() and rule4.empty()):
        if not rule1.empty():
            expand(rule1.pop(), forest, inforest, fullgraph, outsidedegrees, pointers, whichrule, rule1, rule2, rule3, rule4, numFaces)
        elif not rule2.empty():
            expand(rule2.pop(), forest, inforest, fullgraph, outsidedegrees, pointers, whichrule, rule1, rule2, rule3, rule4, numFaces)
        elif not rule3.empty():
            expand(rule3.pop(), forest, inforest, fullgraph, outsidedegrees, pointers, whichrule, rule1, rule2, rule3, rule4, numFaces)
        else:
            expand(rule4.pop(), forest, inforest, fullgraph, outsidedegrees, pointers, whichrule, rule1, rule2, rule3, rule4, numFaces)
    #if our forest does not cover our original graph, we find vertices in our forest that have neighbor(s) outside
    #our forest, and add all its neighbors that are not in the forest. We find such vertices by first finding
    #a vertex not in our forest (if one exists), finding a neighbor of this vertex that is in the forest,
    #and then additing all of this neighbor's neighbors.
    allin=False
    while not allin:
        allin=True
        for key in inforest:
            if not inforest[key]:
                for vertex in fullgraph[key]:
                    if vertex!=-1 and inforest[vertex]:
                        for i in range(numFaces):
                            neighbor=fullgraph[vertex][i]
                            if neighbor!=-1 and not inforest[neighbor]:
                                inforest[neighbor]=True
                                forest[vertex][i]=neighbor
                                if numFaces==6:
                                    forest[neighbor][(i%2+1)%2+i-(i%2)]=vertex
                                else:
                                    forest[neighbor][fullgraph[neighbor].index(vertex)]=vertex
                        break
                if not inforest[key]:
                    allin=False
    #we now turn our forest into a tree, by running bfs (Breadth First Search).
    firstkey=list(forest)[0]
    bfsresult=bfs(forest,firstkey)
    
    #if the maximum distance in the forest is infinite (len(forest) represents infinity)
    if bfsresult['maxdist']==len(forest):
        #the keys that are infinitely far away
        maxkeys=bfsresult['maxkey']
        #while maxkeys isn't empty
        while maxkeys:
            addededge=False
            #for the keys in maxkeys
            for key in maxkeys:
                #for each neighbor in the original graph
                for i in range(numFaces):
                    neighbor=fullgraph[key][i]
                    #checks if the neighbor exists and if it is finitely far away
                    if neighbor!=-1 and bfsresult['distances'][neighbor]<len(forest):
                        #connects the neighbor to the current vertex
                        forest[key][i]=neighbor
                        if numFaces==6:
                            forest[neighbor][(i%2+1)%2+i-(i%2)]=key
                        else:
                            forest[neighbor][fullgraph[neighbor].index(key)]=key
                        #we run bfs again to see if our forest is now a tree
                        bfsresult=bfs(forest,firstkey)
                        #if it is still a forest, we update maxkeys
                        if bfsresult['maxdist']==len(forest):
                            maxkeys=bfsresult['maxkey']
                        #if it is now a tree, we exit
                        else:
                            maxkeys=[]
                        addededge=True
                        break
                if addededge:
                    break
    return forest


# In[67]:

def findSymmetry2(tree, outputfile, numFaces):
    '''
    Algorithm that takes a tree and finds all identical branches within the tree.

    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

    Date: 8/7/2017
        
    Keyword arguments:
    tree -- A directional adjacency list in dictionary form.
    outputfile -- The string path of the file to write the results to
        
    Takes a tree and finds all identical branches within the tree. Does so by first finding the number of
    vertices in every branch in the tree, then comparing the shapes of the branches with the same number
    of vertices. Subsequently, for each edge, eliminates the branch that has fewer matches, then writes
    the information for the top edge of each remaining branch to the output file.
    '''
    #list to hold lists of matched branches
    matches=[]
    #dictionary to hold booleans that state whether one of the two possible trees for an edge has been
    #eliminated.
    edges=dict()
    #dictionary to hold lists of branches where every branch in a given list has the same number of
    #vertices
    bylength=dict()
    count=0
    fraction=0.02
    for vertex in tree:
        #for each neighbor
        for i in range(numFaces):
            #if the neighbor exists
            if tree[vertex][i]!=-1:
                #puts the edge into edges. Note that because we will find each edge twice,
                #we will only add the edge when the parent's id is lower than the child's
                if vertex<tree[vertex][i]:
                    edges[str(vertex)+"-"+str(tree[vertex][i])]=False
                #Uses the explore function to find the number of vertices in the branch
                length=explore2(tree, tree[vertex][i], vertex)
                if length in bylength:
                    bylength[length].append([vertex,i])
                else:
                    bylength[length]=[[vertex,i]]
        count+=1
        if count/fraction>len(tree):
            print(fraction,count)
            fraction+=0.02
    print("Found lengths", len(bylength))
    count=0
    fraction=0.02
    #for each group of branches that have the same length
    for length in bylength:
        numbranches=len(bylength[length])
        #if the length is 1, we know the branches are all matches of each other so we add them to matchees
        #w/o comparing them
        if length==1:
            matches.append([])
            for i in range(numbranches):
                count+=1
                if count>(fraction*len(edges)*2):
                    print(fraction)
                    fraction+=0.02
                matches[len(matches)-1].append([bylength[length][i][0],tree[bylength[length][i][0]][bylength[length][i][1]],[0,1,2,3]])
        else:
            #keeps track of whether we have found match(es) for each branch
            matched=[False]*numbranches
            #for each branch with num vertices = to length
            for i in range(numbranches):
                count+=1
                if not matched[i]:
                    parent1=bylength[length][i][0]
                    direction1=bylength[length][i][1]
                    vertex1=tree[parent1][direction1]
                    #we append this match because it is the first of its type
                    matches.append([[parent1, vertex1, [0]]])
                    matched[i]=True
                    for j in range(i+1,numbranches):
                        if not matched[j]:
                            parent2=bylength[length][j][0]
                            direction2=bylength[length][j][1]
                            vertex2=tree[parent2][direction2]
                            #we compare the two branches using the compare function
                            temp=compare2(tree, vertex1, vertex2, parent1, parent2, direction1, direction2, [True]*4)
                            #if the return from the compare function was [] then there wasn't a match.
                            if temp!=[]:
                                matches[len(matches)-1].append([parent2, vertex2, temp])
                                matched[j]=True
                if count>(fraction*len(edges)*2):
                    print(fraction)
                    fraction+=0.02
    print("Matched")
    #we now sort the matches by the number of branches for each match.
    matchkeys=[]
    matchlengths=[]
    for i in range(len(matches)):
        matchkeys.append(i)
        matchlengths.append(len(matches[i]))
    matchkeys=mergeSort(matchkeys,matchlengths)[0]
    print("Sorted")
    matchpop=[]
    #we get rid of the branch for each edge that has fewer matches
    for key in matchkeys:
        topop=[]
        for i in range(len(matches[key])):
            match=matches[key][i]
            if match[0]<match[1]:
                edgekey=str(match[0])+"-"+str(match[1])
            else:
                edgekey=str(match[1])+"-"+str(match[0])
            if not edges[edgekey]:
                edges[edgekey]=True
                topop.append(i)
        topop=sorted(topop, reverse=True)
        for i in topop:
            matches[key].pop(i)
        if not matches[key]:
            matchpop.append(key)
    print("Found duplicates")
    matchpop=sorted(matchpop, reverse=True)
    #we get rid of any match lists that no longer contain any branches.
    for i in matchpop:
        matches.pop(i)
    print("Eliminated duplicates")
    output=open(outputfile,"w")
    for match in matches:
        for item in match:
            output.write(str(item)+"\t")
        output.write("\n")
    output.close()


# In[69]:

def mergeSort(keys,lengths):
    '''
    Helper function for findSymmetry. Sorts a list of keys by their associated length using merge sort.

    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

    Date: 8/8/2017
        
    Keyword arguments:
    keys -- The keys to be sorted
    lengths -- The lengths associated with the keys.
    '''
    if len(keys)<=1:
        return [keys,lengths]
    leftkeys=[]
    leftlengths=[]
    rightkeys=[]
    rightlengths=[]
    for i in range(len(keys)//2):
        leftkeys.append(keys[i])
        leftlengths.append(lengths[i])
    for i in range(len(keys)//2,len(keys)):
        rightkeys.append(keys[i])
        rightlengths.append(lengths[i])
    left=mergeSort(leftkeys,leftlengths)
    right=mergeSort(rightkeys,rightlengths)
    outkeys=[]
    outlengths=[]
    i=0
    while left[0] and right[0]:
        if left[1][0] <= right[1][0]:
            outkeys.append(left[0].pop(0))
            outlengths.append(left[1].pop(0))
        else:
            outkeys.append(right[0].pop(0))
            outlengths.append(right[1].pop(0))
    while left[0]:
        outkeys.append(left[0].pop(0))
        outlengths.append(left[1].pop(0))
    while right[0]:
        outkeys.append(right[0].pop(0))
        outlengths.append(right[1].pop(0))
    return [outkeys, outlengths]


# In[71]:

def explore2(tree, vertex, parent):
    '''
    Helper function for findSymmetry. Recursively finds the number of vertices in a branch of a tree.

    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

    Date: 8/7/2017
        
    Keyword arguments:
    tree -- A directional adjacency list in dictionary form.
    vertex -- The current vertex to be explored.
    parent -- The vertex's parent
        
    Recursively finds the number of vertices in a branch of a tree. Does so in a DFS manner.
    '''
    #We set height to 1 for the current vertex.
    height=1
    #for each neighbor that exists and isn't the vertex's parent, we call explore again.
    for connection in tree[vertex]:
        if not (connection==-1 or connection==parent):
            height+=explore2(tree, connection, vertex)
    return height


# In[73]:

def compare2(tree, vertex1, vertex2, parent1, parent2, direction1, direction2, rotation, numFaces=6):
    '''
    Helper function for findSymmetry. Recursively compares 2 branches and sees if they have the same shape.

    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

    Date: 8/7/2017
        
    Keyword arguments:
    tree -- The whole tree in dictionary directional adjacency list form.
    vertex1 -- The index of the current vertex from the 1st branch that we are comparing.
    vertex2 -- The index of the current vertex from the 2nd branch that we are comparing.
    parent1 -- The index of vertex1's parent.
    parent2 -- The index of vertex2's parent.
    direction1 -- The direction of the top edge of the 1st branch.
    direction2 -- The direction of the top edge of the 2nd branch.
    rotation -- Boolean list of length 4 that keeps track of whether a rotation of the 2nd branch is possible.
        
    Recursively compares 2 branches in a tree to see if they have the same shape. Does so by recursively: 
      - reorienting every vertex in each branch so that the top edge of each branch is pointing in the 
        positive x direction,
      - comparing the neighbor locations for the vertex from the 1st branch with the neighbor locations of 
        the vertex from the 2nd branch
      - repeating the comparison for the 3 other possible rotations around the x axis for the 2nd vertex.
    '''
    if numFaces==6:
        #reorients vertex1 and vertex2 to ensure that the top edge in their branch is pointing in the positive
        #x direction
        reorient1=reorientVertex2(tree[vertex1], direction1)
        reorient2=reorientVertex2(tree[vertex2], direction2)
        #we first check that the -x and +x neighbors of the 2 vertices are compatible. If they are not, then we can
        #ignore all rotations because we can only rotate around the x axis.
        #for -x and +x:
        for i in range(2):
            #if the neighbor for each vertex exists (-1 means it doesn't) and isn't the vertex's parent (we ignore
            #the parent because we have already checked the parent in the previous recursion)
            if reorient1[i]!=-1 and reorient2[i]!=-1 and reorient1[i]!=parent1 and reorient2[i]!=parent2:
                #we go down a recursion level and compare the two neighbors.
                temp=compare2(tree, reorient1[i], reorient2[i], vertex1, vertex2, direction1, direction2, rotation)
                #if there was a mismatch at a lower recursion level, we return [] to indicate that all rotations
                #have been ruled out.
                if temp==[]:
                    return []
            #if one of the neighbors exists and isn't a parent, but the other doesn't exist, we again return []
            #because all rotations have been ruled out.
            elif not ((reorient1[i]==-1 and reorient2[i]==-1) or (reorient1[i]==parent1 and reorient2[i]==parent2)):
                return []
        #for the 4 possible rotations around the x axis.
        for j in range(4):
            #if rotation j has not yet been ruled out.
            if rotation[j]:
                #for -y, +y, -z, +z
                for i in range(2,6):
                    #if j=0 then we don't need to do any rotating
                    if j==0:
                        #if the neighbor for each vertex exists (-1 means it doesn't) and isn't the vertex's parent
                        if reorient1[i]!=-1 and reorient2[i]!=-1 and reorient1[i]!=parent1 and reorient2[i]!=parent2:
                            #we go down a recursion level and compare the two neighbors, but only for the current
                            #rotation (once the range of possible rotations has been limited at a higher recursion
                            #level, we can't reexpand it at a lower one.)
                            temprot=[False]*4
                            temprot[j]=True
                            temp=compare2(tree, reorient1[i], reorient2[i], vertex1, vertex2, direction1, direction2, temprot)
                            #
                            #if there was a mismatch at a lower recursion level
                            if temp==[]:
                                #we rule out the rotation
                                rotation[j]=False
                                #we break out of the i loop because we know that this j is no longer worth exploring
                                break
                        #if one of the neighbors exists and isn't a parent, but the other doesn't exist
                        elif not ((reorient1[i]==-1 and reorient2[i]==-1) or (reorient1[i]==parent1 and reorient2[i]==parent2)):
                            #we break out of the i loop because we know that this j is no longer worth exploring
                            rotation[j]=False
                            break
                    #if we do need to rotate the 2nd vertex
                    else:
                        #rotates the 2nd vertex according to the +y -> +z -> -y -> -z -> +y... pattern
                        k=i
                        #rotates once, twice, or 3 times depending on j.
                        for l in range(j):
                            if k==2:
                                k=4
                            elif k==3:
                                k=5
                            elif k==4:
                                k=3
                            else:
                                k=2
                        #same as before but this time rotating the 2nd vertex
                        if reorient1[i]!=-1 and reorient2[k]!=-1 and reorient1[i]!=parent1 and reorient2[k]!=parent2:
                            temprot=[False]*4
                            temprot[j]=True
                            temp=compare2(tree, reorient1[i], reorient2[k], vertex1, vertex2, direction1, direction2, temprot)
                            if temp==[]:
                                rotation[j]=False
                                break
                        elif not ((reorient1[i]==-1 and reorient2[k]==-1) or (reorient1[i]==parent1 and reorient2[k]==parent2)):
                            rotation[j]=False
                            break
        #we return all rotations that were found to be valid (or [] if none were)
        output=[]
        for j in range(4):
            if rotation[j]:
                output.append(j)
        return output
    #if this method will be use, need to fix
    elif numFaces==4:
        return []
    else:
        raise ValueError("Unsupported number of faces.")


# In[76]:

def reorientVertex2(vertexlist, direction, numFaces=6):
    '''
    Helper method for findSymmetry. Reorients a vertex so that its branch points in the +x direction.

    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu

    Date: 8/7/2017
        
    Keyword arguments:
    vertexlist -- Ordered adjacency list for the vertex that we want to reorient
    direction -- The actual direction of the branch
        
    Takes a vertex and reorients it so that the top edge of the branch that it is in is pointing in
    the positive x direction. Does so by reordering its adjacency list. Makes sure that the reordering
    process is consistent across coordinates.
    '''
    if numFaces==6:
        #if the original direction of the top edge was -x
        if direction==0:
            output=[vertexlist[1], vertexlist[0], vertexlist[2], vertexlist[3], vertexlist[5], vertexlist[4]]
            return output
        #if the original direction was -y
        elif direction==2:
            output=[vertexlist[3], vertexlist[2], vertexlist[1], vertexlist[0], vertexlist[5], vertexlist[4]]
            return output
        #if the original direction was +y
        elif direction==3:
            output=[vertexlist[2], vertexlist[3], vertexlist[1], vertexlist[0], vertexlist[4], vertexlist[5]]
            return output
        #if the original direction was -z
        elif direction==4:
            output=[vertexlist[4], vertexlist[5], vertexlist[1], vertexlist[0], vertexlist[3], vertexlist[2]]
            return output
        #if the original direction was +z
        elif direction==5:
            output=[vertexlist[5], vertexlist[4], vertexlist[0], vertexlist[1], vertexlist[3], vertexlist[2]]
            return output
        #if the original direction was +x
        return vertexlist.copy()
    elif numFaces==4:
        if direction==1:
            output=[vertexlist[1], vertexlist[2], vertexlist[0], vertexlist[3]]
            return output
        elif direction==2:
            output=[vertexlist[2], vertexlist[3], vertexlist[0], vertexlist[1]]
        elif direction==3:
            output=[vertexlist[3], vertexlist[2], vertexlist[0], vertexlist[2]]
        return vertexlist.copy()
    else:
        raise ValueError("Non-supported number of faces")

if __name__=="__main__":
    # In[77]:
    
    #runs mlst on the small examplegraph
    if 0==1:
        result=mlst(exampleGraph1())
        for key in result:
            print(key,':',result[key])
    
    
    # In[78]:
    
    #run mlst algorithm on the wonderwoman logo test
    if 0==1:
        ingraph=open("/Users/Tristan/Documents/ELM_Project/logo_test_adjacency.txt")
        outfile=open("/Users/Tristan/Documents/ELM_Project/logo_test_MLST.txt","w")
        graph=dict()
        for line in ingraph:
            first=line.find(" ")
            second=first+4
            third=line.find(",",second+1)
            fourth=line.find(",",third+1)
            fifth=line.find(",",fourth+1)
            sixth=line.find(",",fifth+1)
            seventh=line.find(",",sixth+1)
            eighth=line.find("]",seventh+1)
            neighbors=[int(line[second:third]),int(line[third+1:fourth]),int(line[fourth+1:fifth]), int(line[fifth+1:sixth]),int(line[sixth+1:seventh]),int(line[seventh+1:eighth])]
            key=int(line[0:first])
            graph[key]=[]
            for n in neighbors:
                graph[key].append(n)
        result=mlst(graph)
        print(len(result))
        for key in result:
            line=str(key)+" : "+str(result[key])+"\n"
            outfile.write(line)
        ingraph.close()
        outfile.close()
    
    
    # In[79]:
    
    if 0==1:
        ingraph=open("/Users/Tristan/Documents/ELM_Project/coolring2_test_adjacency.txt")
        outfile=open("/Users/Tristan/Documents/ELM_Project/coolring2_test_MLST.txt","w")
        graph=dict()
        for line in ingraph:
            first=line.find(" ")
            second=first+4
            third=line.find(",",second+1)
            fourth=line.find(",",third+1)
            fifth=line.find(",",fourth+1)
            sixth=line.find(",",fifth+1)
            seventh=line.find(",",sixth+1)
            eighth=line.find("]",seventh+1)
            neighbors=[int(line[second:third]),int(line[third+1:fourth]),int(line[fourth+1:fifth]), int(line[fifth+1:sixth]),int(line[sixth+1:seventh]),int(line[seventh+1:eighth])]
            key=int(line[0:first])
            graph[key]=[]
            for n in neighbors:
                graph[key].append(n)
        result=mlst(graph)
        print(len(result))
        for key in result:
            line=str(key)+" : "+str(result[key])+"\n"
            outfile.write(line)
        ingraph.close()
        outfile.close()
    if 0==1:
        ingraph=open("/Users/Tristan/Documents/ELM_Project/test_cylinder_1_test_adjacency.txt")
        outfile=open("/Users/Tristan/Documents/ELM_Project/test_cylinder_1_test_MLST.txt","w")
        graph=dict()
        for line in ingraph:
            first=line.find(" ")
            second=first+4
            third=line.find(",",second+1)
            fourth=line.find(",",third+1)
            fifth=line.find(",",fourth+1)
            sixth=line.find(",",fifth+1)
            seventh=line.find(",",sixth+1)
            eighth=line.find("]",seventh+1)
            neighbors=[int(line[second:third]),int(line[third+1:fourth]),int(line[fourth+1:fifth]), int(line[fifth+1:sixth]),int(line[sixth+1:seventh]),int(line[seventh+1:eighth])]
            key=int(line[0:first])
            graph[key]=[]
            for n in neighbors:
                graph[key].append(n)
        result=mlst(graph)
        print(len(result))
        for key in result:
            line=str(key)+" : "+str(result[key])+"\n"
            outfile.write(line)
        ingraph.close()
        outfile.close()
    if 1==1:
        ingraph=open("/Users/Tristan/Documents/ELM_Project/simpleframe_test_adjacency.txt")
        outfile=open("/Users/Tristan/Documents/ELM_Project/simpleframe_test_MLST.txt","w")
        graph=dict()
        for line in ingraph:
            first=line.find(" ")
            second=first+4
            third=line.find(",",second+1)
            fourth=line.find(",",third+1)
            fifth=line.find(",",fourth+1)
            sixth=line.find(",",fifth+1)
            seventh=line.find(",",sixth+1)
            eighth=line.find("]",seventh+1)
            neighbors=[int(line[second:third]),int(line[third+1:fourth]),int(line[fourth+1:fifth]), int(line[fifth+1:sixth]),int(line[sixth+1:seventh]),int(line[seventh+1:eighth])]
            key=int(line[0:first])
            graph[key]=[]
            for n in neighbors:
                graph[key].append(n)
        result=mlst(graph)
        print(len(result))
        for key in result:
            line=str(key)+" : "+str(result[key])+"\n"
            outfile.write(line)
        ingraph.close()
        outfile.close()
    
    
    # In[80]:
    
    #count the number of leaves for the wonderwoman logo test result
    if 0==1:
        tree=open("/Users/tdaif/Documents/ELM_Project/logo_test_MLST.txt")
        leaves=0
        other=0
        problem=0
        for line in tree:
            count=0
            loc=0
            while line.find("-1",loc)!=-1:
                count+=1
                loc=line.find("-1",loc)+1
            if count<5:
                other+=1
            elif count==5:
                leaves+=1
            else:
                problem+=1
                print(line)
        print(leaves,other,problem)
    
    
    # In[81]:
    
    if 0==1:
        ingraph=open("/Users/tdaif/Documents/ELM_Project/logo_test_adjacency.txt")
        graph=dict()
        for line in ingraph:
            first=line.find(" ")
            second=first+4
            third=line.find(",",second+1)
            fourth=line.find(",",third+1)
            fifth=line.find(",",fourth+1)
            sixth=line.find(",",fifth+1)
            seventh=line.find(",",sixth+1)
            eighth=line.find("]",seventh+1)
            neighbors=[int(line[second:third]),int(line[third+1:fourth]),int(line[fourth+1:fifth]), int(line[fifth+1:sixth]),int(line[sixth+1:seventh]),int(line[seventh+1:eighth])]
            key=int(line[0:first])
            graph[key]=[]
            for n in neighbors:
                graph[key].append(n)
        findSymmetry2(mlst(graph),"/Users/tdaif/Documents/ELM_Project/symmetry_test.txt")
    
    
    # In[82]:
    
    if 0==1:
        m=mlst(exampleGraph1())
        print(m)
        findSymmetry2(m,"/Users/tdaif/Documents/ELM_Project/symmetry_test_example2.txt")
        findSymmetry3(m,"/Users/tdaif/Documents/ELM_Project/symmetry_test_example3.txt")
    
    
    # In[83]:
    
    if 0==1:
        print(mlst(exampleGraphSymmetry()))
    
    
    # In[84]:
    
    if 0==1:
        m=exampleGraphSymmetry()
        findSymmetry2(m,"/Users/tdaif/Documents/ELM_Project/symmetry_test_example_symmetry_2.txt")
        findSymmetry3(m,"/Users/tdaif/Documents/ELM_Project/symmetry_test_example_symmetry_3.txt")
    
    
    # In[85]:
    
    if 0==1:
        symmetry=open("/Users/tdaif/Documents/ELM_Project/symmetry_test.txt")
        counts=[]
        for line in symmetry:
            count=0
            loc=0
            while line.find("\t",loc)!=-1:
                count+=1
                loc=line.find("\t",loc)+1
            counts.append(count)
        plt.hist(counts,log=True)
        plt.show()
    
    
    # In[ ]:



