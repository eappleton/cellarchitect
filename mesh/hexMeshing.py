from Make_Foam_Files import readSTL
import math, random
import numpy as np
import pdb, time
from multiprocessing import Lock, Pipe, Process
import multiprocessing

class HexMeshing():
    def hexMesh(self, triangles, dimensions, blocksize, triangleGrid, tgFactor, numCores):
        self.triangles=triangles
        self.blocksize=blocksize
        self.triangleGrid=triangleGrid
        self.tgFactor=tgFactor
        time1=time.time()
        midx=(dimensions[0]+dimensions[1])/2
        midy=(dimensions[2]+dimensions[3])/2
        midz=(dimensions[4]+dimensions[5])/2
        numx=2*math.ceil((midx-dimensions[0])/blocksize)
        self.numy=2*math.ceil((midy-dimensions[2])/blocksize)
        self.numz=2*math.ceil((midz-dimensions[4])/blocksize)
        self.basey=midy-(self.numy*blocksize/2)
        self.basex=midx-(numx*blocksize/2)
        self.basez=midz-(self.numz*blocksize/2)
        totalBlocks=numx*self.numy*self.numz
        self.pointsPerAxis=math.floor((17*(math.e**(-.015*totalBlocks)))+3)
        grid=dict()
        self.blockIncrement=blocksize/(self.pointsPerAxis+1)
        print("pointsPerAxis", self.pointsPerAxis, totalBlocks)
        index=0
        self.starty=[]
        self.startz=[]
        time2=time.time()
        print("Initialization",time2-time1)
        for y in range(self.numy):
            self.starty.append(y*blocksize+self.basey+self.blockIncrement)
        for z in range(self.numz):
            self.startz.append(z*blocksize+self.basez+self.blockIncrement)
        interval=numx//numCores
        startIndex=0
        assignments=[]
        for i in range(numCores-1):
            assignments.append(range(startIndex,startIndex+interval))
            startIndex+=interval
        assignments.append(range(startIndex,numx))
        processes=[]
        lock=Lock()
        end1,end2=Pipe()
        for i in range(numCores-1):
            processes.append(Process(target=self.parallelizeBlockSearch, args=(lock, end2, assignments[i])))
            processes[i].start()
        for x in assignments[-1]:
            a=x*blocksize+self.basex+self.blockIncrement
            for y in range(self.numy):
                b=self.starty[y]
                for z in range(self.numz):
                    c=self.startz[z]
                    count=0
                    if count==0:
                        for i in range(self.pointsPerAxis):
                            if count==0:
                                for j in range(self.pointsPerAxis):
                                    if count==0:
                                        for k in range(self.pointsPerAxis):
                                            if self.checkPoint(a+i*self.blockIncrement,b+j*self.blockIncrement,c+k*self.blockIncrement):
                                                count+=1
                                                break
                    if count>0:
                        if x not in grid:
                            grid[x]=dict()
                        if y not in grid[x]:
                            grid[x][y]=dict()
                        grid[x][y][z]=index
                        index+=1
        received=1
        while received<numCores:
            count=0
            grid2=end1.recv()
            for x in grid2:
                for y in grid2[x]:
                    for z in grid2[x][y]:
                        grid2[x][y][z]=grid2[x][y][z]+index
                        count+=1
            index+=count
            for x in grid2:
                grid[x]=grid2[x]
            received+=1
        time1=time.time()
        print("Found blocks",time1-time2)
        adjacency=dict()
        for x in sorted(grid.keys()):
            for y in grid[x]:
                for z in grid[x][y]:
                    adjacency[grid[x][y][z]]=[-1,-1,-1,-1,-1,-1]
                    if x-1 in grid and y in grid[x-1] and z in grid[x-1][y]:
                        adjacency[grid[x][y][z]][0]=grid[x-1][y][z]
                        adjacency[grid[x-1][y][z]][1]=grid[x][y][z]
                    if y-1 in grid[x] and z in grid[x][y-1]:
                        adjacency[grid[x][y][z]][2]=grid[x][y-1][z]
                        adjacency[grid[x][y-1][z]][3]=grid[x][y][z]
                    if z-1 in grid[x][y]:
                        adjacency[grid[x][y][z]][4]=grid[x][y][z-1]
                        adjacency[grid[x][y][z-1]][5]=grid[x][y][z]
        time2=time.time()
        print("Made adjacency", time2-time1)
        return adjacency
    
    def parallelizeBlockSearch(self, lock, pipe, assignment):
        grid=dict()
        index=0
        for x in assignment:
            a=x*self.blocksize+self.basex+self.blockIncrement
            for y in range(self.numy):
                b=self.starty[y]
                for z in range(self.numz):
                    c=self.startz[z]
                    count=0
                    if count==0:
                        for i in range(self.pointsPerAxis):
                            if count==0:
                                for j in range(self.pointsPerAxis):
                                    if count==0:
                                        for k in range(self.pointsPerAxis):
                                            if self.checkPoint(a+i*self.blockIncrement,b+j*self.blockIncrement,c+k*self.blockIncrement):
                                                count+=1
                                                break
                    if count>0:
                        if x not in grid:
                            grid[x]=dict()
                        if y not in grid[x]:
                            grid[x][y]=dict()
                        grid[x][y][z]=index
                        index+=1
        lock.acquire()
        pipe.send(grid)
        lock.release()
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
    
    def checkPoint(self, x,y,z):
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
        #pdb.set_trace()
        gridY=math.floor(y*self.tgFactor)
        gridZ=math.floor(z*self.tgFactor)
        if gridY not in self.triangleGrid:
            return False
        if gridZ not in self.triangleGrid[gridY]:
            return False
        tol=self.blocksize*(1e-09)
        np.seterr(all='raise')
        countPos=0
        countNeg=0
        problems=[]
        for t in self.triangleGrid[gridY][gridZ]:
            triangle=self.triangles[t]
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
    
    def makePly(self, adjacency, file_path, name):
        pointList=[]
        faceList=[]
        self.explore(adjacency, pointList, faceList)
        outfile=open(file_path, "w")
        outfile.write("ply\n")
        outfile.write("format ascii 1.0\n")
        outfile.write("obj_info "+name+"\n")
        outfile.write("element vertex "+str(len(pointList))+"\n")
        outfile.write("property float x\n")
        outfile.write("property float y\n")
        outfile.write("property float z\n")
        outfile.write("property uchar red\n")
        outfile.write("property uchar green\n")
        outfile.write("property uchar blue\n")
        outfile.write("element face "+str(len(faceList))+"\n")
        outfile.write("property list uchar int vertex_indices\n")
        outfile.write("end_header\n")
        for line in pointList:
            outfile.write(line)
        for line in faceList:
            outfile.write(line)
        outfile.close()
    
    def add(self, v1,v2):
        return [v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]]
    
    def explore(self, adjacency, pointList, faceList):
        visited=set()
        start=list(adjacency.keys())[0]
        visited.add(start)
        queue=[[start,0.0,0.0,0.0]]
        adjustment=[(-4.0,0.0,0.0), (4.0,0.0,0.0), (0.0,-4.0,0.0), (0.0,4.0,0.0), (0.0,0.0,-4.0), (0.0,0.0,4.0)]
        while queue:
            current=queue.pop()
            color=(random.randint(50,255),random.randint(50,255),random.randint(50,255))
            startPoint=len(pointList)
            self.cubePoints(current[1], current[2], current[3], color, pointList)
            self.face4(0+startPoint,1+startPoint,2+startPoint,3+startPoint,faceList)
            self.face4(5+startPoint,6+startPoint,2+startPoint,1+startPoint,faceList)
            self.face4(4+startPoint,5+startPoint,1+startPoint,0+startPoint,faceList)
            self.face4(7+startPoint,4+startPoint,0+startPoint,3+startPoint,faceList)
            self.face4(6+startPoint,7+startPoint,3+startPoint,2+startPoint,faceList)
            self.face4(7+startPoint,6+startPoint,5+startPoint,4+startPoint,faceList)
            for i,neighbor in enumerate(adjacency[current[0]]):
                if neighbor!=-1 and neighbor not in visited:
                    visited.add(neighbor)
                    queue.append([neighbor]+self.add(current[1:],adjustment[i]))
        if len(visited)!=len(adjacency):
            print("not connected")
                
    
    def cubePoints(self, x,y,z,color, pointlist):
        #top
        self.point(x-1.5,y+1.5,z+1.5, color, pointlist)
        self.point(x+1.5,y+1.5,z+1.5, color, pointlist)
        self.point(x+1.5,y+1.5,z-1.5, color, pointlist)
        self.point(x-1.5,y+1.5,z-1.5, color, pointlist)
        
        #bottom
        self.point(x-1.5,y-1.5,z+1.5, color, pointlist)
        self.point(x+1.5,y-1.5,z+1.5, color, pointlist)
        self.point(x+1.5,y-1.5,z-1.5, color, pointlist)
        self.point(x-1.5,y-1.5,z-1.5, color, pointlist)
    
    def point(self, x,y,z,color,pointlist):
        pointlist.append(str(x)+" "+str(y)+" "+str(z)+" "+str(color[0])+" "+str(color[1])+" "+str(color[2])+" \n")
    
    def face4(self, a,b,c,d,facelist):
        facelist.append("4 "+str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" \n")
    
    def run(self, stl_file, name, maxCells, numCores):
        max_optimization_cycles=10
        triangles=readSTL(stl_file)
        dimensions=triangles[:6]
        triangles=triangles[6]
        triangleGrid, tgFactor=self.triageTriangles(triangles, dimensions)
        blocksize=max(dimensions[1]-dimensions[0], dimensions[3]-dimensions[2], dimensions[5]-dimensions[4])
        blocksize=round(blocksize/2, -math.floor(math.log10(blocksize))+2)
        interval=blocksize/2
        adjacency=self.hexMesh(triangles, dimensions, blocksize, triangleGrid, tgFactor, numCores)
        best=0
        best_adjacency=None
        for i in range(max_optimization_cycles):
            print(blocksize, interval)
            numCells=len(adjacency)*8
            print(numCells)
            if numCells<=maxCells:
                if numCells>.75*maxCells:
                    self.makePly(adjacency, stl_file[:-4]+"_mesh.ply", name)
                    print("Num Cells",numCells)
                    return adjacency
                elif numCells>best:
                    best=numCells
                    best_adjacency=adjacency
                blocksize+=-interval
            else:
                blocksize+=interval
            interval=interval/2
            adjacency=self.hexMesh(triangles, dimensions, blocksize, triangleGrid, tgFactor, numCores)
        self.makePly(best_adjacency, stl_file[:-4]+"_mesh.ply", name)
        print("Num Cells",best)
        return best_adjacency

if __name__=="__main__":
    if 0==1:
        hm=HexMeshing()
        hm.run("examples/sphere1.stl", "sphere1", 2048, 4)
    if 1==1:
        hm=HexMeshing()
        hm.run("examples/simpleframe.STL", "simpleframe", 2048, 4)
    if 0==1:
        hm=HexMeshing()
        hm.run("examples/coolring2.STL", "coolring2", 20000, 4)