import copy, random
import pdb
def freeMovingParts2(mst, fullgraph):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/7/2018
    
    Arguments:
    mst -           minimum spanning tree represented in adjacency list format as a dictionary of lists. Each list contains the same number of elements,
                        with -1 for no connection.
    fullgraph -     adjacency list containing all possible connections, in the same format as the mst
    
    Purpose: Adds edges to a minimum spanning tree to try to make it so that each node has either 1 or >2 connections. The point of this is to
        try to lock in a shape by adding more connections to stabilize it. 
    '''
    mst=copy.deepcopy(mst)
    problemNodes=set()
    #counts number of connections for each node
    for node in mst:
        count=0
        for neighbor in mst[node]:
            if neighbor!=-1:
                count+=1
        if count==2:
            problemNodes.add(node)
    solutions=[]
    solutionLookup=dict()
    #for nodes that have 2 connections, tries to find edges that could be added to graph to give more than 2 connections
    for node in fullgraph:
        if node in problemNodes:
            solutionLookup[node]=[]
            for i in range(len(fullgraph[node])):
                neighbor=fullgraph[node][i]
                if neighbor not in mst[node]:
                    if neighbor in problemNodes:
                        #prioritize edges that fix 2 "problem nodes"
                        solution=[2, node, neighbor, i]
                    else:
                        #de-prioritize edges that only fix current "problem node"
                        solution=[1, node, neighbor, i]
                    solutions.append(solution)
                    solutionLookup[node].append(solution)
    if len(solutions)>0:
        for solution in solutions:
            #de-prioritize edges that block other "solution" edges (since we will add at most 1 edge per "problem node")
            if solution[0]==2:
                count=-1
                for solution2 in solutionLookup[solution[2]]:
                    if solution2[0]>1:
                        count+=1
                #can at most be blocking 3 others
                solution[0]+=-count/4
            #de-prioritize edges that go to leaves (ie nodes with only 1 edge) in the graph.
            elif solution[0]==1:
                count=0
                for x in mst[solution[2]]:
                    if x!=-1:
                        count+=1
                if count==1:
                    solution[0]+=-.25
        #sort in decreasing order
        solutions.sort(reverse=True)
        count=0
        maxCount=len(problemNodes)
        index=0
        while count<maxCount and index<len(solutions):
            solution=solutions[index]
            if solution[1] in problemNodes:
                if solution[2] in problemNodes:
                    mst[solution[1]][solution[3]]=solution[2]
                    for solution2 in solutionLookup[solution[2]]:
                        if solution2[2]==solution[1]:
                            mst[solution[2]][solution2[3]]=solution[1]
                            break
                    problemNodes.remove(solution[1])
                    problemNodes.remove(solution[2])
                    count+=2
                elif solution[0]<=1:
                    mst[solution[1]][solution[3]]=solution[2]
                    for i in range(len(fullgraph[solution[2]])):
                        if fullgraph[solution[2]][i]==solution[1]:
                            mst[solution[2]][i]=solution[1]
                            break
                    problemNodes.remove(solution[1])
                    count+=1
            index+=1
    return mst

#kruskal's leaving out the weighted aspect
def minimumSpanningTree(fullgraph, numFaces, seed=True):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/7/2018
    
    Arguments:
    fullgraph - adjacency list containing all possible connections, as a dictionary of lists. Each list contains the same number of elements,
                    with -1 for no connection.
    numFaces -  the number of faces for each block in the mesh
    seed -      boolean. True if the mst algorithm should use a constant random seed (for reproducability), False is should use generated seed.
    
    Purpose: Generates a minimum spanning tree for a connected graph. Applies Kruskal's algorithm but selects a random edge to add since
        there are no edge weights.
    '''
    if seed:
        random.seed(1)
    sets=dict()
    finalgraph=dict()
    for key in fullgraph:
        sets[key]={key}
        finalgraph[key]=[-1]*numFaces
    edges=[]
    for key in fullgraph:
        for i in range(numFaces):
            neighbor=fullgraph[key][i]
            if neighbor>key:
                edges.append((key, neighbor, i, fullgraph[neighbor].index(key)))
    random.shuffle(edges)
    for edge in edges:
        if edge[0] not in sets[edge[1]]:
            newSet=sets[edge[0]].union(sets[edge[1]])
            for node in newSet:
                sets[node]=newSet
            finalgraph[edge[0]][edge[2]]=edge[1]
            finalgraph[edge[1]][edge[3]]=edge[0]
    return finalgraph