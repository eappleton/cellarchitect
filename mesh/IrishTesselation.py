import math
import pdb, time
from multiprocessing import Process, Lock, Pipe
#import multiprocessing

def createIrish(dimensions, numCores=1):
    '''
    Generates an Irish tesselation as described in "Packing, tiling, and covering with tetrahedra" by JH Conway and S Torquato PNAS 7/11/06 https://doi.org/10.1073/pnas.0601389103
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
    Date: 3/28/19
    
    Arguments:
    dimensions -- the dimensions to make the tesselation
    numCores -- the number of compute cores to use
    
    Returns:
    lis of point ids for each tetrahedron, the coordinates of the points, and the dimensions of the tesselation.
    '''
    #time1=time.time()
    #create all the points used in the tesselation
    points, pointLookup=generatePoints(dimensions)
    #time2=time.time()
    #print("Generated Points",time2-time1)
    #determine all points that are close to each other
    pointAdjacency=createPointAdjacency(points, pointLookup, dimensions)
    #time1=time.time()
    #print("pointAdjacency",time1-time2)
    if numCores>1:
        processes=[]
        startIndex=0
        lock=Lock()
        end1,end2=Pipe()
        #split up points amongst the processes
        interval=len(pointAdjacency)//numCores
        keyAssignments=[]
        allKeys=list(pointAdjacency.keys())
        for i in range(numCores-1):
            keyAssignments.append(allKeys[startIndex:startIndex+interval])
            startIndex+=interval
        keyAssignments.append(allKeys[startIndex:])
        #identify the tetrahedrons
        for i in range(numCores-1):
            processes.append(Process(target=findTets, args=(keyAssignments[i],pointAdjacency, lock, end2)))
            processes[i].start()
            startIndex+=interval
        tets=findTets(keyAssignments[-1], pointAdjacency)
        received=0
        while received<numCores-1:
            received+=1
            tets=tets.union(end1.recv())
        tets=list(tets)
    else:
        #identify the tetrahedrons
        tets=list(findTets(pointAdjacency.keys(),pointAdjacency))
    #time2=time.time()
    #print("Found tets",time2-time1)
    #get rid of any extra points that aren't in any tets, and determine dimensions of tesselation
    dimensions, points=cleanPoints(tets,points)
    #time1=time.time()
    #print("Cleaned points", time1-time2)
    print("Done Irish")
    return tets, points, dimensions

def cleanPoints(tets, points):
    '''
    Removes unused points from the tesselation and determines the tesselation's dimensions
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
    Date: 3/28/19
    
    Arguments:
    tets -- list of the points for each tet
    points -- the coordinates for the points
    
    Returns: the dimensions of the tesselation, as well as the cleaned points
    '''
    #time1=time.time()
    usedPoints=set()
    for tet in tets:
        for point in tet:
            usedPoints.add(point)
    #time2=time.time()
    #print("CP found points",time2-time1)
    tempPoint=usedPoints.pop()
    dimensions=[points[tempPoint][0],points[tempPoint][0],points[tempPoint][1],points[tempPoint][1],points[tempPoint][2],points[tempPoint][2]]
    usedPoints.add(tempPoint)
    #time1=time.time()
    #print("CP initialized dimensions", time1-time2)
    mapping=dict()
    offset=0
    points2=[]
    
    for i, point in enumerate(points):
        if i in usedPoints:
            mapping[i]=i-offset
            points2.append(point)
            if point[0]<dimensions[0]:
                dimensions[0]=point[0]
            elif point[0]>dimensions[1]:
                dimensions[1]=point[0]
            if point[1]<dimensions[2]:
                dimensions[2]=point[1]
            elif point[1]>dimensions[3]:
                dimensions[3]=point[1]
            if point[2]<dimensions[4]:
                dimensions[4]=point[2]
            elif point[2]>dimensions[5]:
                dimensions[5]=point[2]
        else:
            offset+=1
    #time2=time.time()
    #print("CP mapped points", time2-time1)
    #fix indices in the tetrahedrons.
    for i, tet in enumerate(tets):
        tets[i]=frozenset(mapping[x] for x in tet)
    #time1=time.time()
    #print("CP fixed points", time1-time2)
    return dimensions, points2
    

def findTets(keys, pointAdjacency, lock=None, pipe=None):
    '''
    Identifies tetrahedrons in the tesselation
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
    Date: 3/28/19
    
    Arguments:
    keys -- the points to start from
    pointAdjacency -- adjacency list for all the points
    
    Returns: the tetrahedrons (as a set of tuples of point ids)
    '''
    start_time=time.time()
    tets=set()
    for point in keys:
        #explore out from current point, looking for loops of size 4
        found=searchForTets(point, pointAdjacency)
        #print(found)
        for tet in found:
            tets.add(frozenset(tet))
    end_time=time.time()
    print("Done",start_time,end_time, end_time-start_time)
    if lock:
        lock.acquire()
        pipe.send(list(tets))
        lock.release()
    else:
        return tets

def searchForTets(start, adjacency):
    '''
    Identifies tetrahedrons that include the start point
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
    Date: 3/28/19
    
    Arguments:
    start -- the point we are starting from
    pointAdjacency -- adjacency list for all the points
    
    Returns: the tetrahedrons (as a set of tuples of point ids)
    
    Identifies tetrahedrons that include the start point by exploring outward from each point and finding loops of size 4
    '''
    neighbors1=adjacency[start]
    tets=[]
    for neighbor1 in neighbors1:
        neighbors2=adjacency[neighbor1]
        first_intersect=neighbors1.intersection(neighbors2)
        for neighbor2 in first_intersect:
            neighbors3=adjacency[neighbor2]
            second_intersect=neighbors3.intersection(first_intersect)
            for neighbor3 in second_intersect:
                tets.append([start, neighbor1, neighbor2, neighbor3])
    return tets
    
def generatePoints(dimensions):
    '''
    Generates the points for the irish tesselation
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
    Date: 3/28/19
    
    Arguments:
    dimensions -- the specificied dimensions of the tesselation.
    
    Returns: a list of points, and a grid sorting points by location
    '''
    quarter=(-.25,.25)
    half=(-.5,.5)
    points=[]
    pointSet=set()
    count=0
    pointLookup=dict()
    #we buffer the grid (ie add empty cells) in the negative direction in order to make searching the first actual row (for each coordinate) the same as searching any other
    #add x buffer
    pointLookup[dimensions[0]-1]=dict()
    x=dimensions[0]-1
    nextPrint=1000
    #add y buffer and normal y to x buffer
    for y in range(dimensions[2]-2, dimensions[3]+1):
        pointLookup[x][y]=dict()
        #add z buffer and normal z to y buffer and normal y of x buffer
        for z in range(dimensions[4]-2, dimensions[5]+1):
            pointLookup[x][y][z]=[]
    for x in range(dimensions[0], dimensions[1]):
        pointLookup[x]=dict()
        #add y buffer to normal x
        pointLookup[x][dimensions[2]-1]=dict()
        pointLookup[x][dimensions[2]-2]=dict()
        pointLookup[x][dimensions[3]]=dict()
        #add z buffer and normal z to y buffer of normal x
        for y in [dimensions[2]-1, dimensions[2]-2, dimensions[3]]:
            for z in range(dimensions[4]-2, dimensions[5]+1):
                pointLookup[x][y][z]=[]
        for y in range(dimensions[2], dimensions[3]):
            pointLookup[x][y]=dict()
            #add z buffer to normal y/normal x
            pointLookup[x][y][dimensions[4]-1]=[]
            pointLookup[x][y][dimensions[4]-2]=[]
            pointLookup[x][y][dimensions[5]]=[]
            for z in range(dimensions[4], dimensions[5]):
                pointLookup[x][y][z]=[]
                candidates=[]
                points.append((x,y,z))
                pointLookup[x][y][z].append(count)
                count+=1
                #add all the "+" points from the paper
                for q in quarter:
                    for h in half:
                        candidates.append((x, y+q, z+h))
                        candidates.append((x+q, y+h, z))
                        candidates.append((x+h, y, z+q))
                #add all the "1" and "3" points from the paper.
                for h1 in range(2):
                    for h2 in range(2):
                        for h3 in range(2):
                            if math.isclose(((x+y+z+half[h1]+half[h2]+half[h3])*2)%2,1, rel_tol=1e-06, abs_tol=1e-06):
                                candidates.append((x+half[h1],y+half[h2], z+half[h3]))
                for candidate in candidates:
                    if candidate not in pointSet:
                        points.append(candidate)
                        pointSet.add(candidate)
                        pointLookup[math.floor(candidate[0])][math.floor(candidate[1])][math.floor(candidate[2])].append(count)
                        count+=1
        if count>nextPrint:
            print(count, x)
            nextPrint=count+1000
    return points, pointLookup

def distance(point1, point2):
    '''
    Calculates the distance between 2 points
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 3/28/19
    
    Arguments:
    point1 -- the coordinates of the 1st point
    point2 -- the coordinates of the 2nd point
    
    Returns: the distance
    '''
    return ((point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2)**0.5

def createPointAdjacency(points, pointLookup, dimensions):
    '''
    Creates an adjacency list of points in space. Points are considered adjacent if they are within .65 units in 3d space
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 3/28/19
    
    Arguments:
    points -- list of tuples that represent 3d coordinates
    pointLookup -- 3d grid that maps indices of points to their location in 3d space.
    dimensions -- the dimensions of the tesselation
    
    Returns: the adjacency list
    
    Uses the pointLookup grid to efficiently identify points that are adjacent.
    '''
    adjacency=dict()
    for i in range(len(points)):
        adjacency[i]=set()
    #.65 was chosen (somewhat arbitrarily) as follows: Let m=the shortest possible tet edge in the tesselation and let M=the longest possible tet edge
    #in the tesselation. Then M<.65<2*m.
    interaction_length=.65
    #we search through the grid, checking the points in the current cell against other points in the current cell, as well as against the points in half the
    #shell of cells surounding the current cell.
    for x in range(dimensions[0]-1, dimensions[1]-1):
        for y in range(dimensions[2]-1, dimensions[3]):
            for z in range(dimensions[4]-1, dimensions[5]):
                for y2 in range(y-1, y+2):
                    for z2 in range(z-1, z+2):
                        for i in pointLookup[x][y][z]:
                            for j in pointLookup[x+1][y2][z2]:
                                d=distance(points[i],points[j])
                                if d<interaction_length:
                                    adjacency[i].add(j)
                                    adjacency[j].add(i)
                for z2 in range(z-1, z+2):
                    for i in pointLookup[x][y][z]:
                        for j in pointLookup[x][y+1][z2]:
                            d=distance(points[i],points[j])
                            if d<interaction_length:
                                adjacency[i].add(j)
                                adjacency[j].add(i)
                for i in pointLookup[x][y][z]:
                    for j in pointLookup[x][y][z+1]:
                        d=distance(points[i],points[j])
                        if d<interaction_length:
                            adjacency[i].add(j)
                            adjacency[j].add(i)
                #check agains other points in current cell
                for a in range(len(pointLookup[x][y][z])):
                    i=pointLookup[x][y][z][a]
                    for b in range(a+1, len(pointLookup[x][y][z])):
                        j=pointLookup[x][y][z][b]
                        d=distance(points[i],points[j])
                        if d<interaction_length:
                            adjacency[i].add(j)
                            adjacency[j].add(i)
    return adjacency

def gradualSave(coordinates, tetrahedrons):
    byX=dict()
    for i,tet in enumerate(tetrahedrons):
        key=min([coordinates[x][0] for x in tet])
        if key not in byX:
            byX[key]=[i]
        else:
            byX[key].append(i)
    xValues=sorted(list(byX.keys()))
    print(xValues)
    for i in range(len(xValues)):
        currentXs=xValues[:i+1]
        currentTets=[]
        for x in currentXs:
            for tet_id in byX[x]:
                currentTets.append(tet_id)
        currentPoints=set()
        for tet_id in currentTets:
            for point in tetrahedrons[tet_id]:
                currentPoints.add(point)
        mapping=dict()
        offset=0
        currentCoordinates=[]
        for j in range(len(coordinates)):
            if j in currentPoints:
                mapping[j]=j-offset
                currentCoordinates.append(coordinates[j])
            else:
                offset+=1
        tetrahedrons2=[]
        for tet_id in currentTets:
            tet=tetrahedrons[tet_id]
            tetrahedrons2.append([mapping[x] for x in tet])
        saveAsMeshFile(currentCoordinates, tetrahedrons2,"examples/Irish_testing/test"+str(i)+".mesh")
        #plotCgal("examples/Irish_testing/test"+str(i)+".mesh", "IrishMeshTest", False)
        
        

def saveAsMeshFile(coordinates, tetrahedrons, outPath):
    file=open(outPath,"w")
    file.write("MeshVersionFormatted 1\nDimension 3\nVertices\n"+str(len(coordinates))+"\n")
    for coordinate in coordinates:
        file.write(str(coordinate[0])+" "+str(coordinate[1])+" "+str(coordinate[2])+" 1\n")
    file.write("Triangles\n0\nTetrahedra\n"+str(len(tetrahedrons))+"\n")
    for points in tetrahedrons:
        for point in points:
            file.write(str(point+1)+" ")
        file.write("1\n")
    file.write("End\n")
    file.close()
        
if __name__=="__main__":
    if 0==0:
        tets, points, dimensions=createIrish([0,5,0,5,0,5])
        saveAsMeshFile(points, tets, "examples/IrishMeshTest.mesh")
        print(dimensions)
    if 0==1:
        tets, points=createIrish([0,3,0,3,0,3])
        gradualSave(points, tets)
    