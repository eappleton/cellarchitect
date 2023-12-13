import re, math
from numpy import matrix
import pdb
def identifyTets(meshFile):
    mesh=open(meshFile)
    count=-1
    maxCount=-1
    reading1=False
    reading2=False
    triangles=[]
    coordinates=[]
    #interpret Mesh file output from gmsh
    for line in mesh:
        #read triangles
        if reading2:
            count+=1
            if count==0:
                maxCount=int(line[:-1])
            elif count>0:
                spaces=[m.start() for m in re.finditer(" ", line)]
                triangles.append([int(line[spaces[-3]+1:spaces[-2]]),int(line[spaces[-2]+1:spaces[-1]]),int(line[spaces[-1]+1:])])
                if count==maxCount:
                    break
        #read point coordinates
        elif reading1:
            count+=1
            if count==0:
                maxCount=int(line[:-1])
            elif count>0:
                spaces=[m.start() for m in re.finditer(" ", line)]
                coordinates.append([float(line[spaces[-3]+1:spaces[-2]]),float(line[spaces[-2]+1:spaces[-1]]),float(line[spaces[-1]+1:])])
                if count==maxCount:
                    reading1=False
        elif "$Elements" in line:
            reading2=True
            count=-1
        elif "$Nodes" in line:
            reading1=True
    #pdb.set_trace()
    #map all triangle connections
    triangleConnections=[]
    for i in range(maxCount):
        triangleConnections.append(set())
    for i in range(maxCount):
        for j in range(i+1,maxCount):
            shared=0
            for x in triangles[i]:
                if x in triangles[j]:
                    shared+=1
            if shared==2:
                triangleConnections[i].add(j)
                triangleConnections[j].add(i)
    
    #find all tetrahedrons based on triangle connections
    tets=set()
    pdb.set_trace()
    for i in range(maxCount):
        for j in triangleConnections[i]:
            if j>i:
                intersection1=triangleConnections[i].intersection(triangleConnections[j])
                if len(intersection1)==2:
                    tets.add(tuple(sorted(list(intersection1)+[i,j])))
    print(tets)
    tets=list(tets)
    
    #make tetrahedron adjacency list. chirality is preserved
    adjacency=dict()
    for i in range(len(tets)):
        adjacency[i]=[]
    for i in range(len(tets)):
        #find the points in the tetrahedron
        pointIDs=triangles[tets[i][0]]
        for point in triangles[tets[i][1]]:
            if point not in pointIDs:
                pointIDs.append(point)
                break
        #reorder the points so that all tetrahedrons have same chirality
        pointIDs=reorderPoints(coordinates, pointIDs)
        orderedTriangles=reorderTet(tets[i], triangles, pointIDs)
        for triangle in orderedTriangles:
            found=False
            for j in range(len(tets)):
                if i!=j and triangle in tets[j]:
                    adjacency[i].append(j)
                    found=True
                    break
            if not found:
                adjacency[i].append(-1)
    return adjacency
            
        
    
def createVector(point1,point2):
    return (point2[0]-point1[0], point2[1]-point1[1], point2[2]-point2[2])

def cross(vector1, vector2):
    return (vector1[1]*vector2[2]-vector1[2]*vector2[1], vector1[0]*vector2[2]-vector1[2]*vector2[0], vector1[0]*vector2[1]-vector1[1]*vector2[0])
    
def reorderPoints(points, pointIDs):
    v1=createVector(points[pointIDs[1]], points[pointIDs[2]])
    v2=createVector(points[pointIDs[2]], points[pointIDs[3]])
    normal=cross(v1,v2)
    mat1=matrix([[normal[0],v1[0],v2[0]], [normal[1],v1[1],v2[1]], [normal[2],v1[2],v2[2]]])
    mat2=mat1.I
    solution=createVector(points[pointIDs[2]], points[pointIDs[0]])
    mat3=matrix([[solution[0]], [solution[1]], [solution[2]]])
    normal_factor=mat2.dot(mat3)[0]
    if normal_factor>0:
        return pointIDs
    else:
        return [pointIDs[0], pointIDs[2], pointIDs[1], pointIDs[3]]
    
def reorderTet(tet, triangles, pointIDs):
    newTet=[]
    for t in tet:
        if pointIDs[1] in triangles[t] and pointIDs[2] in triangles[t] and pointIDs[3] in triangles[t]:
            newTet.append(t)
    for i in range(2,4):
        for t in tet:
            if pointIDs[0] in triangles[t] and pointIDs[i] in triangles[t] and pointIDs[(i+1)%4+1] in triangles[t]:
                newTet.append(t)
                break
    return newTet
        
if __name__=="__main__":
    if 1==1:
        print(identifyTets("examples/simpleframe_tetmesh0.msh"))