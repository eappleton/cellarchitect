import math, pdb
def makeDevTree(meshType, graph, numFaces, colorBlocks=False):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/6/2018
    
    Arguments:
    meshType:   The type of meshing algorithm used. "Hex" or "Tet"
    graph:      Graph represented by an adjacency list stored as a dictionary of lists. Each list has length the number of sides
                    of the blocks in the mesh, with -1 at a spot in the list where there is no connection.
    numFaces:   int giving the number of faces of the blocks in the mesh.
    
    Purpose: Take a graph of connected blocks in a mesh and find a developmental tree that will lead from one cell to the required population
                of cells to produce the shape
    '''
    #TODO add symmetry identification functionality
    
    ####################################################################################
    ######                          Determine Binders                             ######
    ####################################################################################
    fluorescent_proteins=["eGFP", "EBFP2", "mScarlet", "iRFP670", "tdTomato", "EYFP", "mBanana"]
    fp_index=0
    cells=[]
    homocount=0
    block_reference=dict()
    block_id=0
    heterocount=0
    connections=dict()
    ####################################################################################
    ######                          Hex Mesh                                      ######
    ####################################################################################
    if numFaces==6:
        #cells to be given binders to other blocks on each face.
        facelookup=[(0,5), (2,7), (0,3), (4,7), (3,5), (2,4)]
        #corresponding cells on adjacent block to cells on current block for binders
        facelookup2={(0,5):(2,7), (2,7):(0,5), (0,3):(4,7), (4,7):(0,3), (3,5):(2,4), (2,4):(3,5)}
        for key in graph:
            #under current philosophy, requires 3 hetero protein pairs to make an 8 cell cubic block.
            #1 pair to make a trigonal planar shape with 4 cells, 1 pair to make trigonal planar
            #shape with other 4 cells, and then last pair to attach the 3 tips of the trigonal planar
            #shape to the 3 tips on the other trigonal planar shape. #changed to model from report on 4/8/19
            hetero1=generateProteinName(homocount+heterocount)
            heterocount+=1
            hetero2=generateProteinName(homocount+heterocount)
            heterocount+=1
            hetero3=generateProteinName(homocount+heterocount)
            heterocount+=1
            index=len(cells)
            end_tuple=(block_id,False)
            block_reference[hetero2+"1"]=(-1,dict())
            block_reference[hetero2+"1"][1][hetero3+"2"]=end_tuple
            block_reference[hetero3+"2"]=(-1,dict())
            block_reference[hetero3+"2"][1][hetero2+"1"]=end_tuple
            block_reference[hetero3+"2"][1][hetero1+"2"]=end_tuple
            block_reference[hetero1+"2"]=(-1,dict())
            block_reference[hetero1+"2"][1][hetero3+"2"]=end_tuple
            block_reference[hetero1+"1"]=(-1,dict())
            block_reference[hetero1+"1"][1][hetero3+"1"]=end_tuple
            block_reference[hetero3+"1"]=(-1,dict())
            block_reference[hetero3+"1"][1][hetero1+"1"]=end_tuple
            block_reference[hetero3+"1"][1][hetero2+"2"]=end_tuple
            block_reference[hetero2+"2"]=(-1,dict())
            block_reference[hetero2+"2"][1][hetero3+"1"]=end_tuple
            block_id+=1
            #assigns proteins based on scheme described above
            if colorBlocks:
                #cells.append({hetero1+"1", hetero3+"1", "CycleArrest", fluorescent_proteins[fp_index]})
                #cells.append({hetero1+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                #cells.append({hetero2+"1", hetero3+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                #cells.append({hetero1+"1", hetero3+"1", "CycleArrest", fluorescent_proteins[fp_index]})
                #cells.append({hetero2+"1", hetero3+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                #cells.append({hetero1+"1", hetero3+"1", "CycleArrest", fluorescent_proteins[fp_index]})
                #cells.append({hetero2+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                #cells.append({hetero2+"1", hetero3+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero2+"1", hetero3+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero1+"1", hetero3+"1", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero2+"2", hetero3+"1", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero1+"2", hetero3+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero2+"2", hetero3+"1", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero1+"2", hetero3+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero2+"1", hetero3+"2", "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({hetero1+"1", hetero3+"1", "CycleArrest", fluorescent_proteins[fp_index]})
                fp_index+=1
                if fp_index==len(fluorescent_proteins):
                    fp_index=0
            else:
                #cells.append({hetero1+"1", hetero3+"1", "CycleArrest"})
                #cells.append({hetero1+"2", "CycleArrest"})
                #cells.append({hetero2+"1", hetero3+"2", "CycleArrest"})
                #cells.append({hetero1+"1", hetero3+"1", "CycleArrest"})
                #cells.append({hetero2+"1", hetero3+"2", "CycleArrest"})
                #cells.append({hetero1+"1", hetero3+"1", "CycleArrest"})
                #cells.append({hetero2+"2", "CycleArrest"})
                #cells.append({hetero2+"1", hetero3+"2", "CycleArrest"})
                cells.append({hetero2+"1", hetero3+"2", "CycleArrest"})
                cells.append({hetero1+"1", hetero3+"1", "CycleArrest"})
                cells.append({hetero2+"2", hetero3+"1", "CycleArrest"})
                cells.append({hetero1+"2", hetero3+"2", "CycleArrest"})
                cells.append({hetero2+"2", hetero3+"1", "CycleArrest"})
                cells.append({hetero1+"2", hetero3+"2", "CycleArrest"})
                cells.append({hetero2+"1", hetero3+"2", "CycleArrest"})
                cells.append({hetero1+"1", hetero3+"1", "CycleArrest"})
            #if this block has connections going to it from previously done cells
            if key in connections:
                for connection in connections[key]:
                    #add proteins to correct cells
                    cells[index+connection[1]].add(connection[0])
            #deal with connections to cells that have not been done yet
            for i in range(numFaces):
                neighbor=graph[key][i]
                if neighbor!=-1 and neighbor>key:
                    #current scheme uses 2 proteins to bind hex blocks face to face, one at each end of a diagonal
                    #accross the face to lock in orientation
                    homo1=generateProteinName(homocount+heterocount)
                    homocount+=1
                    homo2=generateProteinName(homocount+heterocount)
                    homocount+=1
                    face=facelookup[i]
                    cells[index+face[0]].add(homo1)
                    cells[index+face[1]].add(homo2)
                    #makes it so the neighboring block can express the other half of the proteins.
                    if neighbor in connections:
                        connections[neighbor].append((homo1, facelookup2[face][0]))
                        connections[neighbor].append((homo2, facelookup2[face][1]))
                    else:
                        connections[neighbor]=[(homo1, facelookup2[face][0]), (homo2, facelookup2[face][1])]
                        
    ####################################################################################
    ######                          Tet Mesh                                      ######
    ####################################################################################
    elif numFaces==4:
        for key in graph:
            #under current philosophy, we use one homobinder to make each tet.
            homo1=generateProteinName(homocount+heterocount)
            homocount+=1
            index=len(cells)
            block_reference[homo1]=(block_id,False)
            block_id+=1
            if colorBlocks:
                cells.append({homo1, "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({homo1, "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({homo1, "CycleArrest", fluorescent_proteins[fp_index]})
                cells.append({homo1, "CycleArrest", fluorescent_proteins[fp_index]})
                fp_index+=1
                if fp_index==len(fluorescent_proteins):
                    fp_index=0
            else:
                cells.append({homo1, "CycleArrest"})
                cells.append({homo1, "CycleArrest"})
                cells.append({homo1, "CycleArrest"})
                cells.append({homo1, "CycleArrest"})
            #adds connections from previously done cells
            if key in connections:
                #adds connections while preserving chirality
                for connection in connections[key]:
                    i=graph[key].index(connection[2])
                    if i==0:
                        cells[index+2].add(connection[1])
                        cells[index+1].add(connection[1])
                        cells[index].add(connection[0])
                    elif i==1:
                        cells[index+1].add(connection[0])
                        cells[index+2].add(connection[1])
                        cells[index+3].add(connection[1])
                    elif i==2:
                        cells[index+2].add(connection[0])
                        cells[index].add(connection[1])
                        cells[index+3].add(connection[1])
                    else:
                        cells[index+3].add(connection[0])
                        cells[index].add(connection[1])
                        cells[index+1].add(connection[1])
            for i in range(numFaces):
                neighbor=graph[key][i]
                if neighbor!=-1 and neighbor>key:
                    #current scheme is to use 2 proteins to bind tet blocks together, one heterobinder and 1 homobinder, with the heterobinders
                    #expressed by 2 of the cells and the homobinder in the 3rd to lock in orientation.
                    hetero1=generateProteinName(homocount+heterocount)
                    heterocount+=1
                    homo1=generateProteinName(homocount+heterocount)
                    homocount+=1
                    #preserves chirality
                    if i==0:
                        cells[index+2].add(hetero1+"1")
                        cells[index+1].add(hetero1+"1")
                        cells[index].add(homo1)
                    elif i==1:
                        cells[index+1].add(homo1)
                        cells[index].add(hetero1+"1")
                        cells[index+3].add(hetero1+"1")
                    elif i==2:
                        cells[index+2].add(homo1)
                        cells[index].add(hetero1+"1")
                        cells[index+3].add(hetero1+"1")
                    else:
                        cells[index+3].add(homo1)
                        cells[index+2].add(hetero1+"1")
                        cells[index+1].add(hetero1+"1")
                    #makes it so the neighboring block can express the other half of the proteins.
                    if neighbor in connections:
                        connections[neighbor].append((homo1, hetero1+"2", key))
                    else:
                        connections[neighbor]=[(homo1, hetero1+"2", key)]
    else:
        raise ValueError("Invalid number of faces")
        
    ####################################################################################
    ######                          Create Tree                                   ######
    ####################################################################################    
        
    #sorts cells based on similarity of the set of proteins they are expressing.
    cells=sortCells(cells)
    #find number of divisions required to get all the cells
    floatDivisions=math.log2(len(cells))
    numDivisions=int(floatDivisions)
    if floatDivisions>numDivisions:
        numDivisions+=1
    #the number of cells that require the full number of divisions (unless the number of cells is a power of 2, some cells will require 1
    #fewer divisions than the others)
    cellsAtLowest=2*len(cells)-2**numDivisions
    #lyaers of cells in the tree
    layers=[[]]
    devTree=[]
    #put all the cells that require the full number of divisions in the bottom layer, and add to the tree
    for i in range(cellsAtLowest):
        layers[0].append(Cell(cells[i],i))
        #-1 for no division
        devTree.append([-1, cells[i]])
    for i in range(numDivisions):
        layer=layers[i]
        layers.append([])
        #pdb.set_trace()
        #we move across current layer 2 cells at a time, determining whether the previous cell needs to divide normally or assymmetrically
        #to produce those 2 cells. We thus move up the tree
        for j in range(0,len(layer),2):
            #finds the proteins expressed in only one of the 2 cells
            dif=layer[j].expressed.symmetric_difference(layer[j+1].expressed)
            if len(dif)==0:
                #if the proteins expressed by the 2 cells are the same, checks that their descendants are the same via the Cell.equals method.
                if len(layer[j].expressed)>0 or layer[j].equals(layer[j+1]):
                    #-2 for standard division
                    devTree.append([-2,layer[j].name, layer[j+1].name])
                else:
                    #-3 for asymmetric division
                    devTree.append([-3,layer[j].name, layer[j+1].name])
            else:
                #-3 for asymmetric division
                devTree.append([-3,layer[j].name, layer[j+1].name])
            layer[j].merge(layer[j+1], len(devTree)-1)
            layers[-1].append(layer[j])
        #deals with case where some cells need to stop dividing before others due to number of cells not being a power of 2. Appends these cells
        #to the 2nd layer of cells.
        if i==0:
            for j in range(cellsAtLowest, len(cells)):
                layers[-1].append(Cell(cells[j],len(devTree)))
                devTree.append([-1, cells[j]])
    return devTree, numDivisions, block_reference

class Cell():
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/6/2018
    
    Purpose: Hold information about what proteins a cell is expressing and what proteins its descendant cells are expressing, while providing
    a means of comparing the cells.
    '''
    def __init__(self, proteins, name):
        '''
        Author: Tristan Daifuku, Harvard Medical School
        Date: 12/6/2018
        
        Arguments:
        proteins:   set of protein names (as strings) that the cell is expressing
        name:       id for the cell
            
        Purpose: Constructor.
        '''
        self.tree=[proteins]
        self.expressed=proteins
        self.name=name
    
    def merge(self, otherCell, name):
        '''
        Author: Tristan Daifuku, Harvard Medical School
        Date: 12/6/2018
        
        Arguments:
        otherCell:  Another Cell object
        name:       id for the merged cell
            
        Purpose: Merges the current cell with another cell. Does not produce a new cell, simply edits the first cell, by combining
        the descendant trees of the 2 cells and clear the expressed proteins set.
        '''
        self.tree=self.tree+otherCell.tree
        self.expressed=set()
        self.name=name
    
    def equals(self, otherCell):
        '''
        Author: Tristan Daifuku, Harvard Medical School
        Date: 12/6/2018
        
        Arguments:
        otherCell:  Another Cell object
            
        Purpose: Tests 2 cells for equality by comparing their descendant trees. If the cells have the same descendants, they are equal.
        '''
        #if the 2 cells have a different number of descendants then they must be non-equal
        if len(self.tree)!=len(otherCell.tree):
            return False
        #indices from the 2nd cell that have been matched to descendants in the first cell.
        taken=set()
        #for each descendant of the first cell
        for proteins in self.tree:
            found=False
            #tries to match to a not already matched descendant of the 2nd cell
            for i in range(len(self.tree)):
                if i not in taken:
                    if otherCell.tree[i]==proteins:
                        taken.add(i)
                        found=True
                        break
            #if there was no match, the cells are not equal
            if not found:
                return False
        return True
    
def sortCells(cells):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/6/2018
    
    Arguments:
    cells: list of sets containing proteins expressed by cells
        
    Purpose: Sorts a list of cells based on the fraction of proteins they share with each other.
    '''
    graph=[]
    #find the ratio of shared proteins to total proteins between each pair of cells.
    for i in range(len(cells)):
        graph.append([0.0]*i)
    for i in range(len(cells)):
        cell=cells[i]
        for j in range(i+1, len(cells)):
            cell2=cells[j]
            ratio=len(cell.intersection(cell2))/len(cell.union(cell2))
            graph[j][i]=ratio
    #I think that finding optimal hamiltonian path is np-complete, so will instead use non-optimal greedy algorithm to sort
    #all the connections by their ratios in decreasing order
    connections=[]
    for i in range(len(graph)):
        for j in range(len(graph[i])):
            connections.append((graph[i][j], i, j))
    connections.sort(reverse=True)
    #make a hamiltonian path by adding edges with highest ratios that do not violate hamiltionian reqs
    newGraph=dict()
    sets=dict()
    for i in range(len(cells)):
        newGraph[i]=[]
        sets[i]={i}
    for connection in connections:
        if len(newGraph[connection[1]])<2 and len(newGraph[connection[2]])<2 and connection[2] not in sets[connection[1]]:
            newGraph[connection[1]].append(connection[2])
            newGraph[connection[2]].append(connection[1])
            newSet=sets[connection[1]].union(sets[connection[2]])
            for node in newSet:
                sets[node]=newSet
            if len(newSet)==len(cells):
                break
    cellOrder=[]
    #find start of hamiltonian path
    cell=0
    prevCell=0
    for i in range(len(cells)):
        if len(newGraph[cell])==1:
            break
        if newGraph[cell][0]==prevCell:
            prevCell=cell
            cell=newGraph[cell][1]
        else:
            prevCell=cell
            cell=newGraph[cell][0]
    #starting at one end of hamiltonian path, follow path and add cells to list to be returned
    prevCell=-1
    for i in range(len(cells)):
        cellOrder.append(cells[cell])
        if i<len(cells)-1:
            if newGraph[cell][0]==prevCell:
                prevCell=cell
                cell=newGraph[cell][1]
            else:
                prevCell=cell
                cell=newGraph[cell][0]
    return cellOrder
                    

def generateProteinName(count):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/6/2018
    
    Arguments:
    count: integer count marking how many previous names have been generated
    
    Purpose: Returns strings, counting upwards in the following format: A-, then AA-AZ, then BA-BZ...
    '''
    alphabet=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    places=0
    output=""
    while count>=26**(places+1):
        places+=1
    while places>-1:
        result=count//(26**places)
        if places>0:
            result+=-1
        output=output+alphabet[result]
        count=count%(26**places)
        places+=-1
    return output