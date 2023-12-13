
# coding: utf-8

# In[7]:

import math, re
import pdb

def fullyConnected(adjacency):
    print(adjacency)
    queue=[list(adjacency.keys())[0]]
    visited=set()
    visited.add(queue[0])
    while queue:
        current=queue.pop()
        visited.add(current)
        for neighbor in adjacency[current]:
            if neighbor!=-1 and (neighbor not in visited):
                visited.add(neighbor)
                queue.append(neighbor)
    return len(visited)==len(adjacency)
    


# In[8]:


def interpretFoam(ownerloc, neighborloc, facesloc, pointsloc, outputloc, edgelength, write=False):
    '''
    Interpret snappyhexmesh output files and create a directional adjacency list.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/27/2017
    
    Keyword arguments:
    ownerloc -- The path to the snappyhexmesh-generated owner file in string format.
    neighborloc -- The path to the snappyhexmesh-generated neighbour file in string format.
    facesloc -- The path to the snappyhexmesh-generated faces file in string format.
    pointsloc -- The path to the snappyhexmesh-generated points file in string format.
    outputloc -- The path to the output file for the adjacency list in string format.
    edgelenth -- The length of the edges of the cubes in the mesh.

    Function that interprets the files created by snappyhexmesh and outputs an adjacency
    list as a text file. The adjacency list contains directionality information. Each list
    for each block in the mesh is of length 6; in order each entry in the list corresponds
    to -x, +x, -y, +y, -z, +z. The mesh must consist of identical cubes. Note: the owner
    file consists of a list of block ids in the mesh indexed by the face ids of the mesh.
    The row number at which a block is listed, offset by the number of rows in the header,
    is the face id. Thus each block may be listed more than once. That being said, each
    face only has one owner - and possibly a neighbor if it is an interior face. The
    neighbor file is like the owner file but lists all the neighbors indexed by the face
    ids. The neighbor is shorter than the owner file because it only contains interior
    faces. The neighbor file matches up with the owner file. The faces file lists the 4
    points that define the corners of each face. the points file lists the coordinates of
    each point in the xyz plane.
    '''
    pdb.set_trace()
    owner=open(ownerloc)
    neighbor=open(neighborloc)
    faces=open(facesloc)
    points=open(pointsloc)
    if write:
        output=open(outputloc,"w")
    #dictionary to store the adjacency list as we create it
    adjacency=dict()
    #dictionary to store the faces of every cube
    facematch=dict()
    started=False
    finished=False
    lineN=''
    numO=0
    numN=0
    count=0
    #loops through the owner file
    for lineO in owner:
        if not started:
            #we skip through the header of the owner file
            if "(" in lineO:
                #print("hi")
                started=True
                lineN=neighbor.readline()
                #we subsequent skip through the header of the neighbor file
                while "(" not in lineN:
                    lineN=neighbor.readline()
        # once we've skipped through the headers, and before we reach the end of the owner file
        elif ")" not in lineO:
            numO=int(lineO)
            #if we haven't yet reached the end of the neighbor file. Note: the neighbor file will always
            #be smaller than the owner file, because the neighbor file only contains interior faces while
            #the owner file contains both interior and exterior faces
            if not finished:
                lineN=neighbor.readline()
                #if we haven't reached the end of the neighbor file
                if ")" not in lineN:
                    numN=int(lineN)
                    #if numO is in adjacency it will be in facematch because we always add to both at once
                    #We add the neighbor (numN) to the adjacency list for the owner, and add the face id
                    #(count) to the list of faces for the owner
                    if numO in adjacency:
                        adjacency[numO].append(numN)
                        facematch[numO].append(count)
                    else:
                        adjacency[numO]=[numN]
                        facematch[numO]=[count]
                    #we do the reverse for the neighbor
                    if numN in adjacency:
                        adjacency[numN].append(numO)
                        facematch[numN].append(count)
                    else:
                        adjacency[numN]=[numO]
                        facematch[numN]=[count]
                else:
                    finished=True
            #if we have reached the end of neighbor, we put -1 for the adjacent cell to indicate that there
            #isn't one.
            if finished:
                if numO in adjacency:
                    adjacency[numO].append(-1)
                    facematch[numO].append(count)
                else:
                    adjacency[numO]=[-1]
                    facematch[numO]=[count]
            #we increment count to keep track of the current face id.
            count+=1
        #when we've reached the end of the owner file
        else:
            break
    #calls the reorderAdjacency function in order to add xyz directionality to the adjacency list 
    coordinates=reorderAdjacency(faces, points, adjacency, facematch, edgelength)
    #writes the adjacency list to the output file
    if write:
        line=""
        count=0
        entries=len(adjacency)
        for key in adjacency:
            count+=1
            line=""
            line=line+str(key)+" : ["
            for connection in adjacency[key]:
                line=line+str(connection)+","
            if count<entries:
                line=line[0:len(line)-1]+"]\n"
            else:
                line=line[0:len(line)-1]+"]"
            output.write(line)
        output.close()
        owner.close()
        neighbor.close()
        faces.close()
        points.close()
    return adjacency, coordinates

# In[20]:
    
def downSize(adjacency, coordinates, gridDimensions, blocksize):
    bigblocksize=2*blocksize
    grid=dict()
    minx=len(adjacency)
    miny=len(adjacency)
    minz=len(adjacency)
    maxx=0
    maxy=0
    maxz=0
    for key in adjacency:
        x=coordinates[key][0]
        y=coordinates[key][1]
        z=coordinates[key][2]
        x=math.floor((x-gridDimensions[0])/bigblocksize)
        y=math.floor((y-gridDimensions[2])/bigblocksize)
        z=math.floor((z-gridDimensions[4])/bigblocksize)
        if x not in grid:
            grid[x]=dict()
            if x<minx:
                minx=x
            elif x>maxx:
                maxx=x
        if y not in grid[x]:
            grid[y]=dict()
            if y<miny:
                miny=y
            elif y>maxy:
                maxy=y
        if z not in grid[x][y]:
            grid[x][y][z]=-1
            if z<minz:
                minz=z
            elif z>maxz:
                maxz=z
    adjacency2=dict()
    coordinates2=dict()
    count=0
    newDimensions=[minx*bigblocksize+gridDimensions[0], maxx*bigblocksize+gridDimensions[1], miny*bigblocksize+gridDimensions[2],
                   maxy*bigblocksize+gridDimensions[3], minz*bigblocksize+gridDimensions[4], maxz*bigblocksize+gridDimensions[5]]
    for x in range(minx, maxx+1):
        for y in range(miny, maxy+1):
            for z in range(minz,maxz+1):
                if x in grid and y in grid[x] and z in grid[x][y]:
                    adjacency2[count]=[-1,-1,-1,-1,-1,-1]
                    grid[x][y][z]=count
                    coordinates2[count]=((x+0.5)*bigblocksize+gridDimensions[0], (y+0.5)*bigblocksize+gridDimensions[2], (z+0.5)*bigblocksize+gridDimensions[4])
                    if x-1 in grid and y in grid[x-1] and z in grid[x-1][y]:
                        adjacency2[count][0]=grid[x-1][y][z]
                        adjacency2[grid[x-1][y][z]][1]=count
                    if y-1 in grid[x] and z in grid[x][y-1]:
                        adjacency2[count][2]=grid[x][y-1][z]
                        adjacency2[grid[x][y-1][z]][3]=count
                    if z-1 in grid[x][y]:
                        adjacency2[count][4]=grid[x][y][z-1]
                        adjacency2[grid[x][y][z-1]][5]=count
                    count+=1
    return adjacency2, coordinates2, newDimensions
        
    
# In[9]:

def faceMean(face, pointList):
    coordinates=[0.0,0.0,0.0]
    for point in face:
        for i in range(3):
            coordinates[i]+=pointList[point][i]
    return (coordinates[0]/4.0, coordinates[1]/4.0, coordinates[2]/4.0)
            
    
def reorderAdjacency(faces,points,adjacency,facematch,edgelength):
    '''
    Helper function for interpretFoam function that adds directionality info to the adjacency list
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/27/2017
    
    Keyword arguments:
    faces -- The snappyhexmesh-generated faces file in string format.
    points -- The snappyhexmesh-generated points file in string format.
    adjacency -- The preliminary adjacency list in dict form created in the first part of interpretFoam
    facematch -- Dictionary containing lists of faces for each block
    edgelenth -- The length of the edges of the cubes in the mesh.

    Function that takes the preliminary adjacency list from the interpretFoam function and correctly
    sorts the list for each block according to whether the adjacent block is in the -x, +x, -y, +y,
    -z, or +z direction (in that order). If there is no adjacent block in that direction, leaves -1.
    '''
    #due to snappyhexmesh doing things like using, for example, both -1.03621e-15 and 0 to represent 0
    #it's necessary to say that values are equal if they fall within a certain distance of each
    #other. 1/10th of the edgelength was arbitrarily chosen to be this distance.
    approximation=edgelength/10.0
    started=False
    pointList=[]
    faceList=[]
    count=0
    entries=0
    prev=""
    entry=[]
    #reads all the points from file and stores them to a list
    for line in points:
        #we skip through the header
        if (not started):
            if "(" in line:
                started=True
                fwdp=[m.start() for m in re.finditer("\(", line)]
                revp=[m.start() for m in re.finditer("\)", line)]
                for i, f in enumerate(fwdp[1:]):
                    r=revp[i]
                    c=line[f+1:r]
                    spaces=[m.start() for m in re.finditer(" ", c)]
                    pointList.append([float(c[:spaces[0]]), float(c[spaces[0]+1:spaces[1]]), float(c[spaces[1]+1:])])
                break
                #once we're at the end of the header, we record the number of entries
                #(provided in the file) and begin storing the info
                #entries=int(prev)
            else:
                prev=line
        #we use the parsePoint method to interpret the info and store as a list of floats
        #which are xyz coordinates
        #elif count < entries:
        #    pointList.append(parsePoint(line))
        #    count+=1
        else:
            break
    started=False
    count=0
    #we loop through the faces file
    faceCoordinates=[]
    for line in faces:
        if (not started):
            #again we go through header and record # entries
            if "(" in line:
                started=True
                entries=int(prev)
            else:
                prev=line
        elif count < entries:
            #parseFace parses the line and returns a list of ints corresponding to the IDs of
            #the points at the corner of each face
            entry=parseFace(line)
            found=False
            faceCoordinates.append(faceMean(entry, pointList))
            #we need to determine whether each face is in the xy, yz, or xz plane to do so, we loop
            #through x, y, and z (i loop) and then loop through each corner point of the face (j loop).
            #As we loop through, we check to see whether all the x, y, or z values for the 4 points are
            #within the approximation window of one of the points. If so, this means that the face
            #points in this direction. 
            for i in range(3):
                if not found:
                    first=pointList[entry[0]][i]
                    found=True
                    for j in range(1,4):
                        if first > (approximation + pointList[entry[j]][i]) or first < (pointList[entry[j]][i]-approximation):
                            found=False
                    #we store the coordinate of the face (0 for x, 1 for y, 2 for z) and the value of this
                    #coordinate in faceList.
                    if found:
                        faceList.append([i, pointList[entry[0]][i]])
            count+=1
        else:
            break
    #for a in faceList:
    #    if len(a)!=2:
    #        print(a)
    #we now reorder the adjacency lists based on the faceList we created
    coordinates=dict()
    for key in adjacency:
        coordinates[key]=[0.0,0.0,0.0]
        temp=[0]*6
        infos=[None]*6
        #we loop through the 6 faces of each block
        for i in range(6):
            #for each face, we pull its info from faceList
            info=faceList[facematch[key][i]].copy()
            c=faceCoordinates[facematch[key][i]]
            for i in range(3):
                coordinates[key][i]+=c[i]
            #if key==2516:
             #   print(info)
            #we append the neighbor block that the face goes to to our info
            info.append(adjacency[key][i])
            #if key==2516:
            #    print(info)
            #sorts the faces based on whether they're x, y, or z faces, and whether they are the positive
            #or negative faces in each coordinate (for instance, there are always 2 x faces. ones will
            #have x coordinate n and one will have x coordinate n+1; we put the n+1 one in the 2nd spot
            #in the info list, and the n one in the 1st spot).
            if infos[2*info[0]] is None:
                infos[2*info[0]]=info
            elif infos[2*info[0]][1]<info[1]:
                infos[2*info[0]+1]=info
            else:
                infos[2*info[0]+1]=infos[2*info[0]]
                infos[2*info[0]]=info
        for i in range(3):
            coordinates[key][i]=coordinates[key][i]/6.0
        adjacency[key]=[]
        #stores the sorted list back into adjacency
        for i in range(6):
            adjacency[key].append(infos[i][2])
    return coordinates


# In[10]:


def parsePoint(line):
    '''
    Helper function for reorderAdjacency function that interprets a line from the points file
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/27/2017
    
    Keyword arguments:
    line -- a string w/ format (num num num) read from the points file

    Takes a line read from the point file, and returns a list of 3 floats that were read from the line
    '''
    first=line.find(" ")
    second=line.find(" ",first+1)
    third=line.find(")",second+1)
    x=line[1:first]
    y=line[first+1:second]
    z=line[second+1:third]
    return [float(x),float(y),float(z)]


# In[11]:


def parseFace(line):
    '''
    Helper function for reorderAdjacency function that interprets a line from the faces file
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/27/2017
    
    Keyword arguments:
    line -- a string w/ format 4(num num num num) read from the faces file

    Takes a line read from the point file, and returns a list of 4 ints that were read from the line
    '''
    first=line.find(" ")
    second=line.find(" ",first+1)
    third=line.find(" ",second+1)
    fourth=line.find(")",third+1)
    a=line[2:first]
    b=line[first+1:second]
    c=line[second+1:third]
    d=line[third+1:fourth]
    return [int(a),int(b),int(c),int(d)]


# In[12]:

if __name__=="__main__":
#runs interpretFoam on the wonderwoman logo
    if 0==1:
        owner="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/logo/1/polyMesh/owner"
        neighbor="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/logo/1/polyMesh/neighbour"
        faces="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/logo/1/polyMesh/faces"
        points="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/logo/1/polyMesh/points"
        output="/Users/Tristan/Documents/ELM_Project/logo_test_adjacency.txt"
        interpretFoam(owner,neighbor,faces,points,output,1)
    if 0==1:
        owner="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/coolring2/1/polyMesh/owner"
        neighbor="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/coolring2/1/polyMesh/neighbour"
        faces="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/coolring2/1/polyMesh/faces"
        points="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/coolring2/1/polyMesh/points"
        output="/Users/Tristan/Documents/ELM_Project/coolring2_test_adjacency.txt"
        interpretFoam(owner,neighbor,faces,points,output,25)
    if 0==1:
        owner="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/test_cylinder_1/1/polyMesh/owner"
        neighbor="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/test_cylinder_1/1/polyMesh/neighbour"
        faces="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/test_cylinder_1/1/polyMesh/faces"
        points="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/test_cylinder_1/1/polyMesh/points"
        output="/Users/Tristan/Documents/ELM_Project/test_cylinder_1_test_adjacency.txt"
        interpretFoam(owner,neighbor,faces,points,output,15)
    if 0==1:
        owner="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/simpleframe/1/polyMesh/owner"
        neighbor="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/simpleframe/1/polyMesh/neighbour"
        faces="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/simpleframe/1/polyMesh/faces"
        points="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/simpleframe/1/polyMesh/points"
        output="/Users/Tristan/Documents/ELM_Project/simpleframe_test_adjacency.txt"
        interpretFoam(owner,neighbor,faces,points,output,2)
    if 0==1:
        owner="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubesphere/1/polyMesh/owner"
        neighbor="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubesphere/1/polyMesh/neighbour"
        faces="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubesphere/1/polyMesh/faces"
        points="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubesphere/1/polyMesh/points"
        output="/Users/Tristan/Documents/ELM_Project/tubesphere_test_adjacency.txt"
        interpretFoam(owner,neighbor,faces,points,output,20)
    if 1==1:
        owner="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubespherebig/1/polyMesh/owner"
        neighbor="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubespherebig/1/polyMesh/neighbour"
        faces="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubespherebig/1/polyMesh/faces"
        points="/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubespherebig/1/polyMesh/points"
        output="/Users/Tristan/Documents/ELM_Project/tubespherebig_test_adjacency.txt"
        interpretFoam(owner,neighbor,faces,points,output,35)

