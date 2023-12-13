from Make_Foam_Files import readSTL
import math, itertools, re, inspect, os, multiprocessing
import numpy as np
import pdb, time
from collections import deque as Deque
from IrishTesselation import createIrish
from TetMeshVisulaization import plotCgal
from array import array
from multiprocessing import Lock, Pipe, Process

class TetMeshing():
    def __init__(self):
        self.x_axis=np.array([1,0,0])
    
    def tetMesh(self, triangles, dimensions, triangleGrid, tgFactor,meshSize, numCores, outpath):
        #generate interior points for each tetrahedron in the original full mesh
        #tetPoints, meshCoordinates, tetrahedrons, span=self.getAllTetPoints("examples/cube1.mesh", dimensions)
        time1=time.time()
        tetPoints, meshCoordinates, tetrahedrons, span=self.getAllTetPoints("examples/IrishMeshTest.mesh", dimensions, meshSize, numCores)
        time2=time.time()
        print("Got Tet Points", time2-time1)
        #arbitrary, used for tolerance in checkPoint.
        blocksize=span/20
        #tetQueue=Queue()
        end1, end2=Pipe()
        lock=Lock()
        #for i in range(len(tetPoints)):
        #    tetQueue.put(i)
        assignments=[]
        interval=len(tetPoints)//numCores
        start_index=0
        for i in range(numCores-1):
            assignments.append(range(start_index,start_index+interval))
            start_index+=interval
        assignments.append(range(start_index,len(tetPoints)))
        processes=[]
        for i in range(numCores-1):
            processes.append(Process(target=self.checkTet, args=(triangles, blocksize, triangleGrid, tgFactor, tetPoints, assignments[i], end2, lock)))
            processes[i].start()
        includedTets=set()
        for tet_key in assignments[-1]:
            points=tetPoints[tet_key]
            for point in points:
                #checks if the point is in the shape. if so, includes the tetrahedron in the final mesh.
                if self.checkPoint(point[0], point[1], point[2], triangles, blocksize, triangleGrid, tgFactor):
                    includedTets.add(tet_key)
                    break
        received=1
        while received<numCores:
            received+=1
            includedTets=includedTets.union(end1.recv())
        time1=time.time()
        print("Found tets to keep",time1-time2)
        #pdb.set_trace()
        #generates an adjacency list for the original full graph
        fullMeshGraph=self.createFullMeshGraph(meshCoordinates, tetrahedrons, meshSize)
        time2=time.time()
        print("Full mesh graph", time2-time1)
        #generates adjacency list for the final graph, preserving chriality (using -1 for no neighbor)
        adjacency=self.createTetGraph(fullMeshGraph, includedTets, tetrahedrons, meshCoordinates)
        time1=time.time()
        print("Found adjacency", time1-time2)
        #adjacency_full=self.createTetGraph(fullMeshGraph, set(fullMeshGraph.keys()), tetrahedrons, meshCoordinates)
        #print(adjacency)
        #self.fillGaps(adjacency, fullMeshGraph, includedTets)
        #write out in .mesh format for visualization
        self.saveAsMeshFile(meshCoordinates, tetrahedrons, adjacency, outpath)
        time2=time.time()
        print("Saved", time2-time1)
        #self.saveAsMeshFile(meshCoordinates, tetrahedrons, adjacency_full, "examples/test_cut_out_full.mesh")
        #re-index mesh components to large gaps in index.
        adjacency=self.fixIndices(adjacency, len(fullMeshGraph))
        time1=time.time()
        print("Reindexed", time1-time2)
        return adjacency
    
    def checkTet(self, triangles, blocksize, triangleGrid, tgFactor, tetPoints, assignment, pipe, lock):
        includedTets=set()
        for tet_key in assignment:
            points=tetPoints[tet_key]
            for point in points:
                #checks if the point is in the shape. if so, includes the tetrahedron in the final mesh.
                if self.checkPoint(point[0], point[1], point[2], triangles, blocksize, triangleGrid, tgFactor):
                    includedTets.add(tet_key)
                    break
        lock.acquire()
        pipe.send(includedTets)
        lock.release()
    
    def saveAsMeshFile(self, coordinates, tetrahedrons, adjacency, outPath):
        validPoints=set()
        for tet in adjacency:
            for point in tetrahedrons[tet]:
                validPoints.add(point)
        coordinates2=[]
        offset=-1
        mapping=dict()
        #reindex points
        for i in range(len(coordinates)):
            if i in validPoints:
                mapping[i]=i-offset
                coordinates2.append(coordinates[i])
            else:
                offset+=1
        file=open(outPath,"w")
        file.write("MeshVersionFormatted 1\nDimension 3\nVertices\n"+str(len(coordinates2))+"\n")
        for coordinate in coordinates2:
            file.write(str(coordinate[0])+" "+str(coordinate[1])+" "+str(coordinate[2])+" 1\n")
        file.write("Triangles\n0\nTetrahedra\n"+str(len(adjacency))+"\n")
        for tet in adjacency:
            points=tetrahedrons[tet]
            for point in points:
                file.write(str(mapping[point])+" ")
            file.write("1\n")
        file.write("End\n")
        file.close()
    
    #if use, should check that works
    def bfs(self, adjacency, start, target=None, maxDistance=None, exactDistance=None):
        queue=Deque()
        queue.append(start)
        visited={start: start}
        distances=dict()
        for x in adjacency:
            distances[x]=-1
        distances[start]=0
        if target:
            for neighbor in adjacency[start]:
                if neighbor==target:
                    return 1
                visited[neighbor]=start
                queue.append(neighbor)
            while len(queue)>0:
                current=queue.popleft()
                distances[current]=distances[visited[current]]+1
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        if neighbor==target:
                            return distances[current]+1
                        visited[neighbor]=current
                        queue.append(neighbor)
            return -1
        elif maxDistance:
            for neighbor in adjacency[start]:
                visited[neighbor]=start
                queue.append(neighbor)
            while len(queue)>0:
                current=queue.popleft()
                if distances[visited[current]]==maxDistance:
                    return distances
                distances[current]=distances[visited[current]]+1
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited[neighbor]=current
                        queue.append(neighbor)
            return distances
        elif exactDistance:
            atCurrent=[]
            for neighbor in adjacency[start]:
                visited[neighbor]=start
                queue.append(neighbor)
            currentDistance=1
            while len(queue)>0:
                current=queue.popleft()
                if distances[visited[current]]==currentDistance:
                    currentDistance+=1
                    if currentDistance>exactDistance:
                        return atCurrent
                    else:
                        atCurrent=[]
                atCurrent.append(current)
                distances[current]=distances[visited[current]]+1
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited[neighbor]=current
                        queue.append(neighbor)
            return []
        else:
            for neighbor in adjacency[start]:
                visited[neighbor]=start
                queue.append(neighbor)
            while len(queue)>0:
                current=queue.popleft()
                distances[current]=distances[visited[current]]+1
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited[neighbor]=current
                        queue.append(neighbor)
            return distances
                    
                
    def createTetGraph(self, fullMeshGraph, includedTets, tetrahedrons, coordinates):
        adjacency=dict()
        for tet in fullMeshGraph:
            if tet in includedTets:
                neighbors=fullMeshGraph[tet]
                #make sure chirality is preserved
                neighbors=self.sortNeighbors(tet, neighbors, tetrahedrons, coordinates)
                adjacency[tet]=[]
                for neighbor in neighbors:
                    if neighbor in includedTets:
                        adjacency[tet].append(neighbor)
                    else:
                        adjacency[tet].append(-1)
        return adjacency
    
    def sortNeighbors(self, current, neighbors, tetrahedrons, coordinates):
        point1=coordinates[tetrahedrons[current][0]]
        point2=coordinates[tetrahedrons[current][1]]
        point3=coordinates[tetrahedrons[current][2]]
        point4=coordinates[tetrahedrons[current][3]]
        names=[tetrahedrons[current][0],tetrahedrons[current][1],tetrahedrons[current][2],tetrahedrons[current][3]]
        # get vector normal to triangle formed by points1-3
        cross=np.cross(np.array([point2[0]-point1[0], point2[1]-point1[1], point2[2]-point1[2]]),
                 np.array([point3[0]-point1[0], point3[1]-point1[1], point3[2]-point1[2]]))
        #determine whether the normal vector is pointing toward or away from point 4. Reorder points 2 and 3 if necessary
        factor=((point1[0]-point4[0])*cross[0]+(point1[1]-point4[1])*cross[1]+(point1[2]-point4[2])*cross[2])/(cross[0]**2+cross[1]**2+cross[2]**2)
        if factor<0:
            temp=names[1]
            names[1]=names[2]
            names[2]=temp
        neighbors2=[-1,-1,-1,-1]
        for neighbor in neighbors:
            points=tetrahedrons[neighbor]
            if names[0] in points:
                if names[1] in points:
                    if names[2] in points:
                        neighbors2[1]=neighbor
                    else:
                        neighbors2[0]=neighbor
                else:
                    neighbors2[2]=neighbor
            else:
                neighbors2[3]=neighbor
        return neighbors2
        
        
    
    def fixIndices(self, adjacency,maxIndex):
        mapping={-1:-1}
        offset=0
        for i in range(maxIndex):
            if i in adjacency:
                mapping[i]=i-offset
            else:
                offset+=1
        adjacency2=dict()
        for tet in adjacency:
            adjacency2[mapping[tet]]=[mapping[x] for x in adjacency[tet]]
        return adjacency2
    
    def fillGaps(self, adjacency, fullMeshGraph, includedTets):
        pdb.set_trace()
        for tet in includedTets:
            distances=self.bfs(adjacency,tet)
            distancesFull=self.bfs(fullMeshGraph, tet, exactDistance=2)
            for tet2 in distancesFull:
                if tet2 in includedTets and (distances[tet2]>4 or distances[tet2]==-1):
                    print("Found")
                    for neighbor in fullMeshGraph[tet]:
                        if neighbor in fullMeshGraph[tet2]:
                            adjacency[neighbor]=[]
                            for neighbor2 in fullMeshGraph[neighbor]:
                                if neighbor2 in adjacency:
                                    adjacency[neighbor2].append(neighbor)
                                    adjacency[neighbor].append(neighbor2)
                            distances=self.bfs(adjacency, tet)
                            break
    
    def createFullMeshGraph(self, coordinates, tetrahedrons, meshSize):
        pyFile=inspect.stack()[-1].filename
        path="Tet_base_meshes/"
        if "tetMeshing" not in pyFile:
            path="mesh/"+path
        files=os.listdir(path)
        targetFile="full_mesh_graph_"+str(meshSize)+".bin"
        if targetFile in files:
            adjacency=dict()
            file=open(path+targetFile, "rb")
            #numCoordinates, numTets, interiorPointsPerTet
            lengthArray=array('i')
            lengthArray.fromfile(file, 1)
            values=array('i')
            values.fromfile(file,lengthArray[0]*5)
            #pdb.set_trace()
            for i in range(0,len(values),5):
                adjacency[values[i]]=[]
                for j in range(i+1, i+5):
                    if values[j]!=-1:
                        adjacency[values[i]].append(values[j]) 
        else:
            pointsToTets=dict()
            for i,tet in enumerate(tetrahedrons):
                for point in tet:
                    if point not in pointsToTets:
                        pointsToTets[point]={i}
                    else:
                        pointsToTets[point].add(i)
            adjacency=dict()
            for i in range(len(tetrahedrons)):
                adjacency[i]=[]
            #pdb.set_trace()
            for i, tet in enumerate(tetrahedrons):
                pointCombos=itertools.combinations(tet, 3)
                for combo in pointCombos:
                    #tets are neighbors if they share 3 points.
                    neighbors=(pointsToTets[combo[0]].intersection(pointsToTets[combo[1]])).intersection(pointsToTets[combo[2]])
                    for neighbor in neighbors:
                        if neighbor>i:
                            adjacency[neighbor].append(i)
                            adjacency[i].append(neighbor)
            forWrite=[]
            for key in adjacency:
                forWrite.append(key)
                count=4
                for neighbor in adjacency[key]:
                    forWrite.append(neighbor)
                    count+=-1
                for i in range(count):
                    forWrite.append(-1)
            writeArray=array('i',[len(adjacency)]+forWrite)
            file=open(path+targetFile, "wb")
            writeArray.tofile(file)
        file.close()
        return adjacency
    
    def readMeshFile(self, meshFile):
        mesh=open(meshFile)
        tetrahedrons=[]
        coordinates=[[]]
        numPoints=0
        while True:
            line=mesh.readline()
            if "Vertices" in line:
                line=mesh.readline()
                numPoints=int(line)
                break
        for i in range(numPoints):
            line=mesh.readline()
            spaces=[m.start() for m in re.finditer(" ", line)]
            coordinates.append([float(line[:spaces[0]]), float(line[spaces[0]+1:spaces[1]]), float(line[spaces[1]+1:spaces[2]])])
        line=mesh.readline()
        if "Triangles" not in line:
            print("Problem")
        numTriangles=int(mesh.readline())
        for i in range(numTriangles):
            mesh.readline()
        line=mesh.readline()
        if "Tetrahedra" not in line:
            print("Problem")
        numTets=int(mesh.readline())
        for i in range(numTets):
            line=mesh.readline()
            spaces=[m.start() for m in re.finditer(" ", line)]
            info=[int(line[:spaces[0]]), int(line[spaces[0]+1:spaces[1]]), int(line[spaces[1]+1:spaces[2]]), int(line[spaces[2]+1:spaces[3]])]
            tetrahedrons.append(info)
        return coordinates, tetrahedrons
    
    def getAllTetPoints(self, meshFile, dimensions, meshSize, numCores):
        span=0.0
        for i in range(3):
            distance=dimensions[2*i+1]-dimensions[2*i]
            if distance>span:
                span=distance
        span=max(dimensions[1]-dimensions[0], dimensions[3]-dimensions[2], dimensions[5]-dimensions[4])
        factor=2
        scale=factor*span/meshSize
        midMesh=(meshSize-1)/2.0
        displacements=((dimensions[1]+dimensions[0])/2-scale*(midMesh-.5),
                       (dimensions[3]+dimensions[2])/2-scale*midMesh,
                       (dimensions[5]+dimensions[4])/2-scale*midMesh)
        pyFile=inspect.stack()[-1].filename
        path="Tet_base_meshes/"
        if "tetMeshing" not in pyFile:
            path="mesh/"+path
        files=os.listdir(path)
        targetFile="base_mesh_"+str(meshSize)+".bin"
        if targetFile in files:
            intValues=array('i')
            file=open(path+targetFile, "rb")
            #numCoordinates, numTets, interiorPointsPerTet
            intValues.fromfile(file, 3)
            coordinateArray=array('d')
            coordinateArray.fromfile(file, 3*intValues[0])
            tetArray=array('i')
            tetArray.fromfile(file, 4*intValues[1])
            interiorPointArray=array('d')
            interiorPointArray.fromfile(file, intValues[1]*intValues[2]*3)
            file.close()
            coordinates=[]
            for i in range(0, len(coordinateArray), 3):
                coordinates.append((coordinateArray[i]*scale+displacements[0],
                                    coordinateArray[i+1]*scale+displacements[1],
                                    coordinateArray[i+2]*scale+displacements[2]))
            tetrahedrons=[]
            for i in range(0, len(tetArray),4):
                tetrahedrons.append([tetArray[i],tetArray[i+1], tetArray[i+2], tetArray[i+3]])
            points=dict()
            key=0
            for i in range(0, len(interiorPointArray),3*intValues[2]):
                points[key]=[]
                for j in range(i, i+3*intValues[2],3):
                    points[key].append([interiorPointArray[j]*scale+displacements[0],
                                       interiorPointArray[j+1]*scale+displacements[1],
                                       interiorPointArray[j+2]*scale+displacements[2]])
                #if key<1000:
                #    print(points[key])
                key+=1
        else:
            if meshSize<5:
                accuracy=3
            elif meshSize<16:
                accuracy=2
            elif meshSize<26:
                accuracy=1
            else:
                accuracy=0
            tetrahedrons, old_coordinates, mesh_dimensions=createIrish([0,meshSize,0,meshSize,0,meshSize], numCores)
            print(mesh_dimensions)
            absTol=1e-7
            points=dict()
            coordinates=[]
            for coordinate in old_coordinates:
                coordinates.append((coordinate[0]*scale+displacements[0], coordinate[1]*scale+displacements[1], coordinate[2]*scale+displacements[2]))
            for i,tet in enumerate(tetrahedrons):
                points[i]=self.generateTetPoints([old_coordinates[tet[0]], old_coordinates[tet[1]], old_coordinates[tet[2]], old_coordinates[tet[3]]], absTol, accuracy)
            #store tesselation as binary file to not have to regenerate each time. Use binary for lossless save of floats.
            file=open(path+targetFile, "wb")
            numInteriorPoints=[1,5,24,63]
            intValues=array('i',[len(coordinates), len(tetrahedrons), numInteriorPoints[accuracy]])
            intValues.tofile(file)
            for coordinate in old_coordinates:
                cArray=array('d',coordinate)
                cArray.tofile(file)
            for tet in tetrahedrons:
                tArray=array('i', tet)
                tArray.tofile(file)
            for key in range(len(tetrahedrons)):
                for point in points[key]:
                    pArray=array('d',point)
                    pArray.tofile(file)
            file.close()
            for key in points:
                for point in points[key]:
                    for i in range(3):
                        point[i]=point[i]*scale+displacements[i]
                #if key<1000:
                #    print(points[key])
            print(len(points),len(coordinates))
            #print(tetrahedrons[:1000])
        return points, coordinates, tetrahedrons, span
        
    def generateTetPoints(self, vertices, absTol, accuracy):
        mid1=[(vertices[0][0]+vertices[1][0])/2,(vertices[0][1]+vertices[1][1])/2,(vertices[0][2]+vertices[1][2])/2]
        mid2=[(vertices[2][0]+vertices[3][0])/2,(vertices[2][1]+vertices[3][1])/2,(vertices[2][2]+vertices[3][2])/2]
        mid3=[(vertices[0][0]+vertices[2][0])/2,(vertices[0][1]+vertices[2][1])/2,(vertices[0][2]+vertices[2][2])/2]
        mid4=[(vertices[1][0]+vertices[3][0])/2,(vertices[1][1]+vertices[3][1])/2,(vertices[1][2]+vertices[3][2])/2]
        #find 2 medians of the tetrahedron
        v1=[mid1[0]-mid2[0], mid1[1]-mid2[1], mid1[2]-mid2[2]]
        v2=[mid3[0]-mid4[0], mid3[1]-mid4[1], mid3[2]-mid4[2]]
        #solve for their intersection, which is the centroid
        mat=np.array([[v1[0],-v2[0]],[v1[1],-v2[1]]])
        if not math.isclose(np.linalg.det(mat),0.0,rel_tol=1e-06, abs_tol=absTol):
            solution=np.array([mid4[0]-mid2[0],mid4[1]-mid2[1]])
        else:
            mat=np.array([[v1[0],-v2[0]],[v1[2],-v2[2]]])
            if not math.isclose(np.linalg.det(mat),0.0,rel_tol=1e-06, abs_tol=absTol):
                solution=np.array([mid4[0]-mid2[0],mid4[2]-mid2[2]])
            else:
                mat=np.array([[v1[1],-v2[1]],[v1[2],-v2[2]]])
                solution=np.array([mid4[1]-mid2[1],mid4[2]-mid2[2]])
        values=np.linalg.solve(mat, solution)
        centroid=[mid2[0]+v1[0]*values[0],mid2[1]+v1[1]*values[0],mid2[2]+v1[2]*values[0]]
        #24 points
        if accuracy==2:
            points=[]
            mid5=[(vertices[0][0]+vertices[3][0])/2,(vertices[0][1]+vertices[3][1])/2,(vertices[0][2]+vertices[3][2])/2]
            mid6=[(vertices[2][0]+vertices[1][0])/2,(vertices[2][1]+vertices[1][1])/2,(vertices[2][2]+vertices[1][2])/2]
            #determine check points based on centroid.
            for exteriorPoint in vertices+[mid1,mid2,mid3,mid4, mid5, mid6]:
                points.append([(exteriorPoint[0]-centroid[0])*.5+centroid[0],(exteriorPoint[1]-centroid[1])*.5+centroid[1], (exteriorPoint[2]-centroid[2])*.5+centroid[2]])
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            triangle_centroids=[]
            triangle_centroids.append([(2/3)*mid1[0]+vertices[2][0]/3, (2/3)*mid1[1]+vertices[2][1]/3, (2/3)*mid1[2]+vertices[2][2]/3])
            triangle_centroids.append([(2/3)*mid2[0]+vertices[1][0]/3, (2/3)*mid2[1]+vertices[1][1]/3, (2/3)*mid2[2]+vertices[1][2]/3])
            triangle_centroids.append([(2/3)*mid4[0]+vertices[0][0]/3, (2/3)*mid4[1]+vertices[0][1]/3, (2/3)*mid4[2]+vertices[0][2]/3])
            triangle_centroids.append([(2/3)*mid3[0]+vertices[3][0]/3, (2/3)*mid3[1]+vertices[3][1]/3, (2/3)*mid3[2]+vertices[3][2]/3])
            for exteriorPoint in triangle_centroids:
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            return points
        #1 point
        elif accuracy==0:
            return [centroid]
        #5 points
        elif accuracy==1:
            points=[]
            for exteriorPoint in vertices:
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            points.append(centroid)
            return points
        #63 points
        elif accuracy==3:
            points=[centroid]
            mid5=[(vertices[0][0]+vertices[3][0])/2,(vertices[0][1]+vertices[3][1])/2,(vertices[0][2]+vertices[3][2])/2]
            mid6=[(vertices[2][0]+vertices[1][0])/2,(vertices[2][1]+vertices[1][1])/2,(vertices[2][2]+vertices[1][2])/2]
            #determine check points based on centroid.
            for exteriorPoint in vertices+[mid1,mid2,mid3,mid4, mid5, mid6]:
                points.append([(exteriorPoint[0]-centroid[0])*.3+centroid[0],(exteriorPoint[1]-centroid[1])*.3+centroid[1], (exteriorPoint[2]-centroid[2])*.3+centroid[2]])
                points.append([(exteriorPoint[0]-centroid[0])*.6+centroid[0],(exteriorPoint[1]-centroid[1])*.6+centroid[1], (exteriorPoint[2]-centroid[2])*.6+centroid[2]])
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            triangle_centroids=[]
            triangle_centroids.append([(2/3)*mid1[0]+vertices[2][0]/3, (2/3)*mid1[1]+vertices[2][1]/3, (2/3)*mid1[2]+vertices[2][2]/3])
            triangle_centroids.append([(2/3)*mid2[0]+vertices[1][0]/3, (2/3)*mid2[1]+vertices[1][1]/3, (2/3)*mid2[2]+vertices[1][2]/3])
            triangle_centroids.append([(2/3)*mid4[0]+vertices[0][0]/3, (2/3)*mid4[1]+vertices[0][1]/3, (2/3)*mid4[2]+vertices[0][2]/3])
            triangle_centroids.append([(2/3)*mid3[0]+vertices[3][0]/3, (2/3)*mid3[1]+vertices[3][1]/3, (2/3)*mid3[2]+vertices[3][2]/3])
            triangle_centroids.append([(1/3)*mid3[0]+vertices[3][0]*(2/3), (1/3)*mid3[1]+vertices[3][1]*(2/3), (1/3)*mid3[2]+vertices[3][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid6[0]+vertices[3][0]*(2/3), (1/3)*mid6[1]+vertices[3][1]*(2/3), (1/3)*mid6[2]+vertices[3][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid1[0]+vertices[3][0]*(2/3), (1/3)*mid1[1]+vertices[3][1]*(2/3), (1/3)*mid1[2]+vertices[3][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid6[0]+vertices[0][0]*(2/3), (1/3)*mid6[1]+vertices[0][1]*(2/3), (1/3)*mid6[2]+vertices[0][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid4[0]+vertices[0][0]*(2/3), (1/3)*mid4[1]+vertices[0][1]*(2/3), (1/3)*mid4[2]+vertices[0][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid2[0]+vertices[0][0]*(2/3), (1/3)*mid2[1]+vertices[0][1]*(2/3), (1/3)*mid2[2]+vertices[0][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid3[0]+vertices[1][0]*(2/3), (1/3)*mid3[1]+vertices[1][1]*(2/3), (1/3)*mid3[2]+vertices[1][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid2[0]+vertices[1][0]*(2/3), (1/3)*mid2[1]+vertices[1][1]*(2/3), (1/3)*mid2[2]+vertices[1][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid5[0]+vertices[1][0]*(2/3), (1/3)*mid5[1]+vertices[1][1]*(2/3), (1/3)*mid5[2]+vertices[1][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid1[0]+vertices[2][0]*(2/3), (1/3)*mid1[1]+vertices[2][1]*(2/3), (1/3)*mid1[2]+vertices[2][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid4[0]+vertices[2][0]*(2/3), (1/3)*mid4[1]+vertices[2][1]*(2/3), (1/3)*mid4[2]+vertices[2][2]*(2/3)])
            triangle_centroids.append([(1/3)*mid5[0]+vertices[2][0]*(2/3), (1/3)*mid5[1]+vertices[2][1]*(2/3), (1/3)*mid5[2]+vertices[2][2]*(2/3)])
            for exteriorPoint in triangle_centroids:
                points.append([(exteriorPoint[0]-centroid[0])*.5+centroid[0],(exteriorPoint[1]-centroid[1])*.5+centroid[1], (exteriorPoint[2]-centroid[2])*.5+centroid[2]])
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            return points
        return []
                
    #sort triangles into a y-z grid. Allows for quick lookup of relevant triangles when checking if a point is in the shape.
    def triageTriangles(self, triangles, dimensions):
        factor=1000.0/max(dimensions[3]-dimensions[2], dimensions[5]-dimensions[4])
        grid=dict()
        for i,triangle in enumerate(triangles):
            ymin=math.floor(factor*min(triangle[0][1], triangle[1][1], triangle[2][1]))
            ymax=math.ceil(factor*max(triangle[0][1], triangle[1][1], triangle[2][1]))
            zmin=math.floor(factor*min(triangle[0][2], triangle[1][2], triangle[2][2]))
            zmax=math.ceil(factor*max(triangle[0][2], triangle[1][2], triangle[2][2]))
            for y in range(ymin, ymax+1):
                if y not in grid:
                    grid[y]=dict()
                for z in range(zmin, zmax+1):
                    if z not in grid[y]:
                        grid[y][z]=[]
                    grid[y][z].append(i)
        return grid, factor
    
    def checkPoint(self, x,y,z,triangles,blocksize, triangleGrid, tgFactor):
        '''
        Helper function for findInteriorPoint. Checks if a point is inside the stl shape.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 8/1/2017
        
        Keyword arguments:
        x -- The x coordinate of the point of interest as a float.
        y -- The y coordinate of the point of interest as a float.
        z -- The z coordinate of the point of interest as a float.
        triangles -- A 3D list of all the triangles described in the stl.
        blocksize -- The edge length of the blocks in the blockMesh.
    
        Checks if a point is inside the stl shape by counting how many surfaces in the shape
        the positive x vector <+infinity,0,0> from the point (x,y,z) crosses. 
        Attemps to account for floating point error.
        '''
        #if 13<x<22 and 45<y<50 and 45<z<50:
        #    pdb.set_trace()
        #pdb.set_trace()
        gridY=math.floor(y*tgFactor)
        gridZ=math.floor(z*tgFactor)
        if gridY not in triangleGrid:
            return False
        if gridZ not in triangleGrid[gridY]:
            return False
        tol=blocksize*(1e-09)
        np.seterr(all='raise')
        countPos=0
        countNeg=0
        problems=[]
        for t in triangleGrid[gridY][gridZ]:
            triangle=triangles[t]
            tooclose=False
            for problem in problems:
                for vertex in problem:
                    if vertex in triangle:
                        tooclose=True
                    else:
                        tooclose=False
                        break
                if tooclose:
                    break
            #if the positive x vector doesn't go through a vertex of the triangle and if we haven't
            #counted it
            if not tooclose:
                #pdb.set_trace()
                #finds a vector that is normal to the plane that the triangle ABC is in using the cross product of
                #vector A-B and vector A-C
                line1=[triangle[0][0]-triangle[1][0],triangle[0][1]-triangle[1][1],triangle[0][2]-triangle[1][2]]
                line2=[triangle[0][0]-triangle[2][0],triangle[0][1]-triangle[2][1],triangle[0][2]-triangle[2][2]]
                normal=np.cross(line1,line2)
                if np.dot((1,0,0),normal)!=0:
                    w=(x-triangle[0][0],y-triangle[0][1],z-triangle[0][2])
                    x_intersect=x+(-normal[0]*w[0]-normal[1]*w[1]-normal[2]*w[2])/normal[0]
                    xclose=math.isclose(x_intersect, x, rel_tol=1e-06, abs_tol=tol)
                    #we can't safely tell which side of the triangle the point is on.
                    if xclose:
                        return False
                    isPositive=x_intersect>x
                    #if the vector intersects the plane, we see if the point of intersection is in the triangle.
                    #Barycentric technique.
                    #We need to solve the equation <A to P>=u*<A to B>+v*<A to C>, where <> indicate a
                    #vector, P is our point, and ABC are the vertices of our triangle. Note that this is
                    #equivalent to solving a system of equations of the following type:
                    #x_P-x_A=u*(x_B-x_A)+v*(x_C-x_A)
                    #y_P-y_A=u*(y_B-y_A)+v*(y_C-y_A)
                    #Thus:
                    #line1=[x_B-x_A, x_C-x_A]
                    line1=[triangle[1][0]-triangle[0][0], triangle[2][0]-triangle[0][0]]
                    #line2=[y_B-y_A, y_C-y_A]
                    line2=[triangle[1][1]-triangle[0][1], triangle[2][1]-triangle[0][1]]
                    #line3=[z_B-z_A, z_C-z_A]
                    line3=[triangle[1][2]-triangle[0][2], triangle[2][2]-triangle[0][2]]
                    #Note that we need line3 in case the matrix created by line1 and line2
                    #is singular. It is impossible for all 3 combinations of lines to give
                    #singular matrices because ABC is a triangle. We thus find a combination
                    #that gives a non-singular matrix by checking the determinants.
                    if not math.isclose((line1[0]*line2[1]), (line1[1]*line2[0]), rel_tol=1e-06, abs_tol=tol):
                        solution=[x_intersect-triangle[0][0],y-triangle[0][1]]
                        intest=np.linalg.solve(np.array([line1,line2]),np.array(solution))
                    elif not math.isclose((line1[0]*line3[1]), (line1[1]*line3[0]), rel_tol=1e-06, abs_tol=tol):
                        solution=[x_intersect-triangle[0][0],z-triangle[0][2]]
                        intest=np.linalg.solve(np.array([line1,line3]),np.array(solution))
                    else:
                        solution=[y-triangle[0][1],z-triangle[0][2]]
                        intest=np.linalg.solve(np.array([line2,line3]),np.array(solution))
                    #if the point is in the triangle
                    if (intest[0]>0 or math.isclose(intest[0],0, rel_tol=1e-06, abs_tol=tol)) and (intest[1]>0 or math.isclose(intest[1],0, rel_tol=1e-06, abs_tol=tol)) and ((intest[0]+intest[1])<1 or math.isclose(intest[0]+intest[1],1, rel_tol=1e-06, abs_tol=tol)):
                        #otherwise, we deal with cases where the point is very close to the edge of a triangle
                        if math.isclose(intest[0],0, rel_tol=1e-06, abs_tol=tol) or math.isclose(intest[1],0, rel_tol=1e-06, abs_tol=tol) or math.isclose(intest[0]+intest[1],1, rel_tol=1e-06, abs_tol=tol):
                            #deals with axis hitting corner of triangle
                            if math.isclose(intest[0],0, rel_tol=1e-06, abs_tol=tol) and math.isclose(intest[1],0, rel_tol=1e-06, abs_tol=tol):
                                problems.append([triangle[0]])
                            elif math.isclose(intest[0],0, rel_tol=1e-06, abs_tol=tol) and math.isclose(intest[1],1, rel_tol=1e-06, abs_tol=tol):
                                problems.append([triangle[2]])
                            elif math.isclose(intest[0],1, rel_tol=1e-06, abs_tol=tol) and math.isclose(intest[1],0, rel_tol=1e-06, abs_tol=tol):
                                problems.append([triangle[1]])
                            #deals with axis hitting edge of triangle
                            elif math.isclose(intest[0],0, rel_tol=1e-06, abs_tol=tol):
                                problems.append([triangle[0], triangle[2]])
                            elif math.isclose(intest[1],0, rel_tol=1e-06, abs_tol=tol):
                                problems.append([triangle[0], triangle[1]])
                            else:
                                problems.append([triangle[1], triangle[2]])
                        #no matter what we add 1 to the count
                        if isPositive:
                            countPos+=1
                        else:
                            countNeg+=1
        #if the count is odd then the point is inside the shape.
        if countPos%2 == 1 and countNeg%2==1:
            return True
        return False
    
    def run(self, stl_file, name, maxCells, numCores):
        outpath=stl_file[:-4]+"_tet_mesh.mesh"
        triangles=readSTL(stl_file)
        dimensions=triangles[:6]
        triangles=triangles[6]
        triangleGrid, tgFactor=self.triageTriangles(triangles, dimensions)
        meshSize=2
        adjacency=self.tetMesh(triangles, dimensions, triangleGrid, tgFactor, meshSize, numCores, outpath)
        last=adjacency
        while len(adjacency)*4<maxCells:
            meshSize+=1
            last=adjacency
            adjacency=self.tetMesh(triangles, dimensions, triangleGrid, tgFactor, meshSize, numCores, outpath)
        adjacency=last
        print("Num cells", len(adjacency)*4)
        plotCgal(outpath, name, False)
        return adjacency

if __name__=="__main__":
    multiprocessing.set_start_method('fork')
    tm=TetMeshing()
    #print(tm.generateTetPoints([[5,0,0],[0,5,0],[-1,-1,0],[1,1,5]],1e-7))
    #tm.run("examples/simpleframe.STL", "simple_frame_cutout", 256, 4)
    print(tm.run("examples/sphere1.stl", "sphere1", 512, 4))
    #plotCgal("examples/test_cut_out_full.mesh", "test_cut_out_full", False)