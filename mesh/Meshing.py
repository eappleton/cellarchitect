import math, random, itertools, inspect, os, multiprocessing
import numpy as np
import pdb, time
from multiprocessing import Lock, Pipe, Process
from IrishTesselation import createIrish
from array import array

class Meshing():
    '''
    Master class for Meshing, both Tet Meshing and Hex Meshing classes inherit from this class. Contains methods used by both.
    '''
    def readSTL(self, stlPath):
        '''
        Reads an STL file and returns its extreme xyz coordinates and its triangles.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 7/31/2017
        
        Keyword arguments:
        stlPath -- String containing the location on the disc of the STL file.
    
        Reads either a binary or ASCII STL file and returns the extreme xyz coordinates and a list containing all of
        the triangles that made up the STL file. Also, if the STL file was binary, returns the path to a copy of the
        file where the 80 bit header was written over.
        '''
        stl=open(stlPath)
        xmin=0
        xmax=0
        ymin=0
        ymax=0
        zmin=0
        zmax=0
        text=True
        firstentry=True
        triangles=[]
        path=None
        try:
            for line in stl:
                text=True
        except UnicodeDecodeError:
            text=False
        #if the file is ASCII
        if text:
            #reopens the file as a text file
            stl.close()
            stl2=open(stlPath)
            count=0
            triangle=[]
            for line in stl2:
                try:
                    if "vertex" in line:
                        first=0
                        for i in range(len(line)):
                            if line[i].isdigit() or line[i]=="-":
                                first=i
                                break
                        second=line.find(" ",first+1)
                        third=line.find(" ",second+1)
                        x=float(line[first:second])
                        y=float(line[second+1:third])
                        z=float(line[third+1:])
                        triangle.append([x,y,z])
                        if count == 2:
                            triangles.append(triangle)
                            triangle=[]
                            count=0
                        else:
                            count+=1
                        if firstentry:
                            xmin=x
                            xmax=x
                            ymin=y
                            ymax=y
                            zmin=z
                            zmax=z
                            firstentry=False
                        else:
                            if xmin > x:
                                xmin = x
                            elif xmax < x:
                                xmax = x
                            if ymin > y:
                                ymin = y
                            elif ymax < y:
                                ymax = y
                            if zmin > z:
                                zmin=z
                            elif zmax < z:
                                zmax=z
                except ValueError:
                    pdb.set_trace()
            stl2.close()
        #if the file is binary
        else:
            stl.close()
            stl2=open(stlPath,"rb")
            #reads off the 80 bit header
            stl2.read(80)
            #reads the number of triangles and uses interpretInt to convert it to an int.
            numtri=self.interpretInt(stl2.read(4))
            #loops through all the triangles
            for i in range(numtri):
                stl2.read(12)
                triangle=[]
                for j in range(3):
                    x=self.interpretFloat(stl2.read(4))
                    y=self.interpretFloat(stl2.read(4))
                    z=self.interpretFloat(stl2.read(4))
                    triangle.append([x,y,z])
                    if firstentry:
                        xmin=x
                        xmax=x
                        ymin=y
                        ymax=y
                        zmin=z
                        zmax=z
                        firstentry=False
                    else:
                        if xmin > x:
                            xmin = x
                        elif xmax < x:
                            xmax = x
                        if ymin > y:
                            ymin = y
                        elif ymax < y:
                            ymax = y
                        if zmin > z:
                            zmin=z
                        elif zmax < z:
                            zmax=z
                triangles.append(triangle)
                stl2.read(2)
            stl2.close()
            #path=fixSTL(stlPath)
        return [xmin,xmax,ymin,ymax,zmin,zmax,triangles,path]

    def interpretFloat(self, fourbytes):
        '''
        Helper function for readSTL that takes a 4 byte float read from an STL file and returns a python float.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 7/31/2017
        
        Keyword arguments:
        fourbytes -- bytearray of length 4 read from an STL file that represents a float
    
        Takes a bytearray of length 4 read from an STL file and returns a python float. The sign bit
        is the 25th bit (1st bit of last byte), the exponent is bits 26-32,17 (bits 2-7 of late byte
        followed by the 1st bit of the 3rd byte), and the significand is bits 18-24,9-16,1-8 (bits
        2-7 of the 3rd byte, the 2nd byte, then the 1st byte).
        '''
        ints=[]
        for byte in fourbytes:
            ints.append(byte)
        #checks if the sign bit is 1. if so, the float is negative
        if ints[3]>=128:
            sign=-1
            #subtracts out the sign bit and shifts the remaining bits in the byte left by one (since
            #they are bits 1-7 of the exponent)
            exp=(ints[3]-128)*2
        else:
            sign=1
            #shifts the remaining bits in the byte left by one (since they are bits 1-7 of the exponent)
            exp=ints[3]*2
        #determines the value of the 8th bit of the exponent
        if ints[2]>=128:
            exp+=1
            ints[2]+=-128
        signif=(ints[2]*65536)+(ints[1]*256)+ints[0]
        return sign*(1+(signif/8388608))*(2**(exp-127))


    def interpretInt(self, fourbytes):
        '''
        Helper function for readSTL that takes a 4 byte positive int read from an STL file and returns a python int.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 7/31/2017
        
        Keyword arguments:
        fourbytes -- bytearray of length 4 read from an STL file that represents a positive integer
    
        Takes a bytearray of length 4 read from an STL file and returns a python int. The bytes are in reverse order.
        '''
        return fourbytes[0]+(fourbytes[1]*256)+(fourbytes[2]*65536)+(fourbytes[3]*16777216)

    #sort triangles into a y-z grid. Allows for quick lookup of relevant triangles when checking if a point is in the shape.
    def triageTriangles(self):
        '''
        Sort triangles from the input STL file into a y-z grid. Allows for quick lookup of relevant triangles when checking if a point is in the shape.
        
        Author: Tristan Daifuku
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/2019
        '''
        #We make a 1000x1000 grid (if the grid isn't square, then the larger dimension is 1000). 
        #Factor maps the coordinates in the shape to the indices of the grid.
        factor=1000.0/max(self.dimensions[3]-self.dimensions[2], self.dimensions[5]-self.dimensions[4])
        grid=dict()
        for i,triangle in enumerate(self.triangles):
            #puts the triangle's id into every cell in the grid that overlaps the triangle
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
        self.triangleGrid=grid
        self.tgFactor=factor
    
    def checkPoint(self, x,y,z):
        '''
        Checks if a point is inside the stl shape.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Keyword arguments:
        x -- The x coordinate of the point of interest as a float.
        y -- The y coordinate of the point of interest as a float.
        z -- The z coordinate of the point of interest as a float.
    
        Checks if a point is inside the stl shape by counting how many surfaces in the shape
        the vectors <+infinity,0,0> and <-infinity,0,0> from the point (x,y,z) crosses.
        If both numbers are odd, the point is determined to be inside the shape. Uses the triangle
        grid generated by triageTriangles to do this efficiently: points are only compared to triangles
        that are in the same general y-z region of space. Attemps to account for floating point error
        (primarily introduced by the ascii version of stl files)
        '''
        #pdb.set_trace()
        #grid indices of the point
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
            #accounts for situations where a point is so close to one of the sides of the triangle, that it might be confused as being
            #in the neighboring triangle. Should occur very rarely, if at all.
            tooclose=False
            for problem in problems:
                for vertex in problem:
                    if vertex in triangle:
                        tooclose=True
                    else:
                        tooclose=False
                        break
                if tooclose:
                    #print("Too Damn Close")
                    break
            if not tooclose:
                #pdb.set_trace()
                #finds a vector that is normal to the plane that the triangle ABC is in using the cross product of
                #vector A-B and vector A-C
                line1=[triangle[0][0]-triangle[1][0],triangle[0][1]-triangle[1][1],triangle[0][2]-triangle[1][2]]
                line2=[triangle[0][0]-triangle[2][0],triangle[0][1]-triangle[2][1],triangle[0][2]-triangle[2][2]]
                normal=np.cross(line1,line2)
                #checks that the plane of the triangle is not parallel to the x axis (if it is then there is no way for our x vector to cross it)
                if np.dot((1,0,0),normal)!=0:
                    #determines x value at which the x vector crosses the triangles plane. Uses method described on this site: 
                    #http://geomalgorithms.com/a05-_intersect-1.html Does so by finding w, which is a vector from point A
                    #on the triangle to the point we are testing. Then finds the scalar k by which we can multiply <1,0,0> by to have
                    #w+<k,0,0> dot n=0$ with n the normal line found previously (essentially, the vector from A to the point of intersection
                    #must be perpendicular to the normal). The x_intersection thus occurs at x_original-k (we already know the y and z
                    #coordinates of the intersect since our vector only changes in the x direction).
                    w=(x-triangle[0][0],y-triangle[0][1],z-triangle[0][2])
                    x_intersect=x+(-normal[0]*w[0]-normal[1]*w[1]-normal[2]*w[2])/normal[0]
                    #Checks if the point is too close to the plane to be able to tell which side of the triangle the point is on (because of fp error)
                    xclose=math.isclose(x_intersect, x, rel_tol=1e-06, abs_tol=tol)
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
        #if both counts are odd then the point is inside the shape. We check both counts for increased accuracy/consistency.
        if countPos%2 == 1 and countNeg%2==1:
            return True
        return False
    
    def run(self, stl_file, name, maxCells, numCores):
        '''
        Sets up meshing run.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        stl_file -- the path to the stl input file.
        name -- some name associated with the shape.
        maxCells -- integer maximum number of biological cells that can be used in the shape.
        numCores -- the number of cores on the current compute system to use.
    
        Provides framework for setting up meshing runs. Reads in the stl file, and creates the triangle grid used to check points.
        '''
        self.numCores=numCores
        self.triangles=self.readSTL(stl_file)
        self.dimensions=self.triangles[:6]
        self.triangles=self.triangles[6]
        self.triageTriangles()
        
    def point(self, x,y,z,color,pointlist):
        '''
        Makes a point in the ply format. Here because it is used by both Tet and Hex meshing
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        x -- the x coordinate of the point
        y -- the y coordinate of the point
        z -- the z coordinate of the point
        color -- the color of the point as a tuple of ints
        pointlist -- the list of point codes for the ply file
        '''
        pointlist.append(str(x)+" "+str(y)+" "+str(z)+" "+str(color[0])+" "+str(color[1])+" "+str(color[2])+" \n")

    def isFullyConnected(self, adjacency):
        '''
        Checks if a mesh is fully connected (ie there are no free floating elements)
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 4/1/19
        
        Arguments:
        adjacency -- adjacency list stored as a dictionary, with -1 meaning no neighbor on a given face of the block.
        '''
        start=list(adjacency.keys())[0]
        visited={start}
        stack=[start]
        while len(stack)>0:
            current=stack.pop()
            for neighbor in adjacency[current]:
                if neighbor!=-1 and neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        if len(visited)!=len(adjacency):
            raise NotConnectedError
            
class NotConnectedError(Exception):
    '''
    Exception that is raised if a mesh is not fully connected.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 4/1/19
    '''
    pass
        

class HexMeshing(Meshing):
    '''
    Class that produces a mesh of a shape using cube shaped blocks. Does so by trimming blocks not in the shape. Extends Meshing class.
    '''
    def run(self, stl_file, name, maxCells, numCores):
        '''
        Generates a hex mesh based on a stl file.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        stl_file -- the path to the stl input file.
        name -- some name associated with the shape.
        maxCells -- integer maximum number of biological cells that can be used in the shape.
        numCores -- the number of cores on the current compute system to use.
        
        Returns: Dictionary containing an adjacency list of the mesh
        
        Runs a hex mesh based on a stl file. does so by trimming blocks not in the shape. can do multiple runs to try to optimize the number of cells
        used.
        '''
        #at 10 cycles, blocksize change will be 2.4e-4*size of shape so more iterations seems unlikely to accomplish much
        max_optimization_cycles=10
        #calls Meshing.run() to begin setup
        super().run(stl_file, name, maxCells, numCores)
        self.midx=(self.dimensions[0]+self.dimensions[1])/2
        self.midy=(self.dimensions[2]+self.dimensions[3])/2
        self.midz=(self.dimensions[4]+self.dimensions[5])/2
        self.blocksize=max(self.dimensions[1]-self.dimensions[0], self.dimensions[3]-self.dimensions[2], self.dimensions[5]-self.dimensions[4])
        self.blocksize=round(self.blocksize/2, -math.floor(math.log10(self.blocksize))+2)
        interval=self.blocksize/2
        adjacency=self.hexMesh()
        best=0
        best_adjacency=None
        #binary search for best blocksize. if can't find optimal, uses closest to optimal
        for i in range(max_optimization_cycles):
            print(self.blocksize, interval)
            numCells=len(adjacency)*8
            print(numCells)
            if numCells<=maxCells:
                #if we are close to our target, we end
                if numCells>.75*maxCells:
                    self.isFullyConnected(adjacency)
                    self.makePly(adjacency, stl_file[:-4]+"_mesh.ply", name)
                    print("Num Cells",numCells)
                    return adjacency
                elif numCells>best:
                    best=numCells
                    best_adjacency=adjacency
                self.blocksize+=-interval
            else:
                self.blocksize+=interval
            interval=interval/2
            adjacency=self.hexMesh()
        self.isFullyConnected(best_adjacency)
        self.makePly(best_adjacency, stl_file[:-4]+"_mesh.ply", name)
        print("Num Cells",best)
        return best_adjacency

    def hexMesh(self):
        '''
        Runs the meshing algorithm using the current blocksize and the shape read in from the stl file.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Returns: Dictionary containing an adjacency list of the mesh.
        
        Runs the meshing algorithm using the current blocksize and the shape read in from the stl file. Determines if blocks should be kept by
        checking points within each block against the target shape (checks to see if point is in shape). Parallelizes this process.
        '''
        time1=time.time()
        #the number of blocks in each direction
        numx=2*math.ceil((self.midx-self.dimensions[0])/self.blocksize)
        self.numy=2*math.ceil((self.midy-self.dimensions[2])/self.blocksize)
        self.numz=2*math.ceil((self.midz-self.dimensions[4])/self.blocksize)
        totalBlocks=numx*self.numy*self.numz
        self.pointsPerAxis=math.floor((17*(math.e**(-.015*totalBlocks)))+3)
        grid=dict()
        self.blockIncrement=self.blocksize/(self.pointsPerAxis+1)
        print("pointsPerAxis", self.pointsPerAxis, totalBlocks)
        index=0
        self.starty=[]
        self.startz=[]
        time2=time.time()
        print("Initialization",time2-time1)
        #store y and z values so we don't have to recalculate each time
        for y in range(-self.numy//2, self.numy//2):
            self.starty.append(y*self.blocksize+self.midy+self.blockIncrement)
        for z in range(-self.numz//2, self.numz//2):
            self.startz.append(z*self.blocksize+self.midz+self.blockIncrement)
        #splits up x values amongst the cores
        interval=numx//self.numCores
        startIndex=-numx//2
        assignments=[]
        for i in range(self.numCores-1):
            assignments.append(range(startIndex,startIndex+interval))
            startIndex+=interval
        assignments.append(range(startIndex,numx//2))
        #starts parallel processes to quickly go through all x values
        processes=[]
        lock=Lock()
        end1,end2=Pipe()
        for i in range(self.numCores-1):
            processes.append(Process(target=self.parallelizeBlockSearch, args=(lock, end2, assignments[i])))
            processes[i].start()
        #do the same actions in current process as in others
        for x in assignments[-1]:
            a=x*self.blocksize+self.midx+self.blockIncrement
            for y in range(self.numy):
                b=self.starty[y]
                for z in range(self.numz):
                    c=self.startz[z]
                    count=0
                    if count==0:
                        #goes through all interior points in block until it finds one in the shape
                        for i in range(self.pointsPerAxis):
                            if count==0:
                                for j in range(self.pointsPerAxis):
                                    if count==0:
                                        for k in range(self.pointsPerAxis):
                                            if self.checkPoint(a+i*self.blockIncrement,b+j*self.blockIncrement,c+k*self.blockIncrement):
                                                count+=1
                                                break
                    #if one was found in the shape
                    if count>0:
                        #we create a 3d grid of included blocks.
                        if x not in grid:
                            grid[x]=dict()
                        if y not in grid[x]:
                            grid[x][y]=dict()
                        grid[x][y][z]=index
                        index+=1
        #gets results from the other processes
        received=1
        while received<self.numCores:
            count=0
            grid2=end1.recv()
            for x in grid2:
                for y in grid2[x]:
                    for z in grid2[x][y]:
                        grid2[x][y][z]=grid2[x][y][z]+index
                        count+=1
            index+=count
            #merges the grid from the other process into the grid from this process.
            for x in grid2:
                grid[x]=grid2[x]
            received+=1
        time1=time.time()
        print("Found blocks",time1-time2)
        adjacency=dict()
        #goes through the created grid and creates an adjacency list
        for i in range(index):
            adjacency[i]=[-1,-1,-1,-1,-1,-1]
        for x in grid:
            for y in grid[x]:
                for z in grid[x][y]:
                    #checks neighbor in -x direction
                    if x-1 in grid and y in grid[x-1] and z in grid[x-1][y]:
                        adjacency[grid[x][y][z]][0]=grid[x-1][y][z]
                        adjacency[grid[x-1][y][z]][1]=grid[x][y][z]
                    #checks neighbor in -y direction
                    if y-1 in grid[x] and z in grid[x][y-1]:
                        adjacency[grid[x][y][z]][2]=grid[x][y-1][z]
                        adjacency[grid[x][y-1][z]][3]=grid[x][y][z]
                    #checks neighbor in -z direction
                    if z-1 in grid[x][y]:
                        adjacency[grid[x][y][z]][4]=grid[x][y][z-1]
                        adjacency[grid[x][y][z-1]][5]=grid[x][y][z]
        time2=time.time()
        print("Made adjacency", time2-time1)
        return adjacency

    def parallelizeBlockSearch(self, lock, pipe, assignment):
        '''
        Checks all blocks in a certain x range to see if they should be kept. Used to parallelize the process of checking blocks.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        lock -- Lock for communicating with the main process
        pipe -- one end of the pipe for communicating with the main process
        assignment -- the x values that the process needs to analyze.
        '''
        #same as in hexMesh method
        grid=dict()
        index=0
        for x in assignment:
            a=x*self.blocksize+self.midx+self.blockIncrement
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
        #send results back to main process
        lock.acquire()
        pipe.send(grid)
        lock.release()

    def makePly(self, adjacency, file_path, name):
        '''
        Creates a ply file based on the adjacency list found by the meshing process. Used to visualize the mesh.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        adjacency -- the adjacency list generated by the meshing process.
        file_path -- the file to write to
        name -- name to give the shape
        
        Spaces out the blocks so that mesh is more visible.
        '''
        pointList=[]
        faceList=[]
        #determine where blocks should be in space relative to each other, and create the points and faces for the ply file
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
        '''
        Adds 2 vectors stored as lists.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        v1 -- the first vector
        v2 -- the second vector
        '''
        return [v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]]
    
    def explore(self, adjacency, pointList, faceList):
        '''
        explores the adjacency list of the mesh, starting from an arbitrary block, and creates the ply components relative to that starting block
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        adjacency -- the adjacency list generated by the meshing process.
        pointList -- list of point codes for the ply file
        faceList -- list of face codes for the ply file
        '''
        visited=set()
        start=list(adjacency.keys())[0]
        visited.add(start)
        #stack implemented with a list
        queue=[[start,0.0,0.0,0.0]]
        #puts 1/3 side length of spacing between blocks
        adjustment=[(-4.0,0.0,0.0), (4.0,0.0,0.0), (0.0,-4.0,0.0), (0.0,4.0,0.0), (0.0,0.0,-4.0), (0.0,0.0,4.0)]
        while queue:
            current=queue.pop()
            color=(random.randint(50,255),random.randint(50,255),random.randint(50,255))
            startPoint=len(pointList)
            #make points for the cube
            self.cubePoints(current[1], current[2], current[3], color, pointList)
            #make the faces of the cube
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
        '''
        Makes the points of a cube in the ply format. uses an edge length of 3
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        x -- the x coordinate of the center of the cube
        y -- the y coordinate of the center of the cube
        z -- the z coordinate of the center of the cube
        color -- the color of the cube as a tuple of ints
        pointlist -- the list of point codes for the ply file
        '''
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
    
    def face4(self, a,b,c,d,facelist):
        '''
        Makes a 4 cornered face in the ply format using provided point ids
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        a -- the first point id
        b -- the second point id
        c -- the third point id
        d -- the 4th point id
        facelist -- the list of face codes for the ply file
        '''
        facelist.append("4 "+str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" \n")
        
class TetMeshing(Meshing):
    '''
    Class that produces a mesh of a shape using tetrahedron shaped blocks. Does so by trimming blocks not in the shape. Extends Meshing class.
    Uses an inital full mesh of an irish tesselation as described in in "Packing, tiling, and covering with tetrahedra" by JH Conway and S 
    Torquato PNAS 7/11/06 https://doi.org/10.1073/pnas.0601389103
    '''
    
    def run(self, stl_file, name, maxCells, numCores):
        '''
        Generates a tet mesh based on a stl file.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        stl_file -- the path to the stl input file.
        name -- some name associated with the shape.
        maxCells -- integer maximum number of biological cells that can be used in the shape.
        numCores -- the number of cores on the current compute system to use.
        
        Returns: Dictionary containing an adjacency list of the mesh
        
        Creates a tet mesh based on a stl file. does so by trimming blocks in an irish tesselation not in the shape. can do multiple runs to try 
        to optimize the number of cells used.
        '''
        outpath=stl_file[:-4]+"_tet_mesh.ply"
        #calls Meshing.run to begin setup and make the triangle grid based on the stl shape
        super().run(stl_file, name, maxCells, numCores)
        self.span=max(self.dimensions[1]-self.dimensions[0], self.dimensions[3]-self.dimensions[2], self.dimensions[5]-self.dimensions[4])
        #2 is the smallest viable size for the irish tesselation. We start at the lowest and work up, because analyzing the lower sizes takes
        #far less time
        self.meshSize=2
        self.factor=3
        adjacency, includedTets=self.tetMesh()
        last=adjacency
        last_included=includedTets
        last_tets=self.tetrahedrons
        last_coords=self.meshCoordinates
        if len(adjacency)*4>maxCells:
            while len(adjacency)*4>maxCells:
                self.factor+=1
                adjacency, includedTets=self.tetMesh()
                last=adjacency
                last_tets=self.tetrahedrons
                last_coords=self.meshCoordinates
                last_included=includedTets
            interval=(self.factor-3)/4
            self.factor=(self.factor+3)/2
            max_optimization_cycles=10
            best=len(adjacency)*4
            for i in range(max_optimization_cycles):
                adjacency, includedTets=self.tetMesh()
                interval=interval/2
                numCells=len(adjacency)*4
                if numCells<=maxCells:
                    #if we are close to our target, we end
                    if numCells>.75*maxCells:
                        last=adjacency
                        last_tets=self.tetrahedrons
                        last_coords=self.meshCoordinates
                        last_included=includedTets
                        break
                    elif numCells>best:
                        last=adjacency
                        last_tets=self.tetrahedrons
                        last_coords=self.meshCoordinates
                        last_included=includedTets
                        best=numCells
                    self.factor+=-interval
                else:
                    self.factor+=interval
        else:
            #we keep increasing the mesh size (thus increasing the number of blocks) until we hit the cell cap
            self.factor=2
            while len(adjacency)*4<maxCells:
                self.meshSize+=1
                last=adjacency
                last_included=includedTets
                last_tets=self.tetrahedrons
                last_coords=self.meshCoordinates
                adjacency, includedTets=self.tetMesh()
        adjacency=last
        self.isFullyConnected(adjacency)
        print("Num cells", len(adjacency)*4)
        #make the ply file for visualization
        self.makePly(outpath, name, last_tets, last_coords, last_included, False)
        return adjacency
    
    def tetMesh(self):
        '''
        Runs the tet meshing algorithm on the shape using the current mesh size.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Returns: Dictionary containing an adjacency list of the mesh, a set of all the tetrahedron ids that are in the mesh.
        '''
        #generate interior points for each tetrahedron in the original full mesh
        time1=time.time()
        self.getAllTetPoints()
        time2=time.time()
        print("Got Tet Points", time2-time1)
        #arbitrary, used for tolerance in checkPoint.
        self.blocksize=self.span/20
        #for communicating between processes
        end1, end2=Pipe()
        lock=Lock()
        #Divy up tetrahedrons in original full tesselation amongst the processes.
        assignments=[]
        interval=len(self.tetPoints)//self.numCores
        start_index=0
        for i in range(self.numCores-1):
            assignments.append(range(start_index,start_index+interval))
            start_index+=interval
        assignments.append(range(start_index,len(self.tetPoints)))
        processes=[]
        #check whether each tetrahedron is in the shape
        for i in range(self.numCores-1):
            processes.append(Process(target=self.checkTet, args=(assignments[i], end2, lock)))
            processes[i].start()
        includedTets=set()
        for tet_key in assignments[-1]:
            points=self.tetPoints[tet_key]
            for point in points:
                #checks if the point is in the shape. if so, includes the tetrahedron in the final mesh.
                if self.checkPoint(point[0], point[1], point[2]):
                    includedTets.add(tet_key)
                    break
        received=1
        #get results from other processes
        while received<self.numCores:
            received+=1
            includedTets=includedTets.union(end1.recv())
        time1=time.time()
        print("Found tets to keep",time1-time2)
        #pdb.set_trace()
        #generates an adjacency list for the original full graph
        fullMeshGraph=self.createFullMeshGraph()
        time2=time.time()
        print("Full mesh graph", time2-time1)
        #generates adjacency list for the final graph, preserving chriality (using -1 for no neighbor)
        adjacency=self.createTetGraph(fullMeshGraph, includedTets)
        self.reorder(adjacency)
        time1=time.time()
        print("Found adjacency", time1-time2)
        #adjacency_full=self.createTetGraph(fullMeshGraph, set(fullMeshGraph.keys()))
        #self.makePly("examples/full_tet_graph.ply", "full_tet_graph", self.tetrahedrons, self.meshCoordinates, set(fullMeshGraph.keys()), False)
        #re-index mesh components to large gaps in index.
        #reindexes adjacency so that there are no gaps in the tet ids.
        adjacency=self.fixIndices(adjacency, len(fullMeshGraph))
        time2=time.time()
        print("Reindexed", time2-time1)
        return adjacency, includedTets
    
    def checkTet(self, assignment, pipe, lock):
        '''
        Checks whether tetrahedrons are in the shape. Used for multiprocessing
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        lock -- Lock for communicating with the main process
        pipe -- one end of the pipe for communicating with the main process
        assignment -- the x values that the process needs to analyze.
        '''
        #does same as in tetMeshing()
        includedTets=set()
        for tet_key in assignment:
            points=self.tetPoints[tet_key]
            for point in points:
                #checks if the point is in the shape. if so, includes the tetrahedron in the final mesh.
                if self.checkPoint(point[0], point[1], point[2]):
                    includedTets.add(tet_key)
                    break
        lock.acquire()
        pipe.send(includedTets)
        lock.release()

    def createTetGraph(self, fullMeshGraph, includedTets):
        '''
        Creates the adjacency list of the mesh based on the included tetrahedrons.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        fullMeshGraph -- adjacency list (w/o positional info) for the full tesselation
        includedTets -- set of tet ids that are included in the final mesh
        
        Returns: adjacency list in dictionary format of the final mesh, with -1 at positions with no neighbor
        '''
        adjacency=dict()
        for tet in fullMeshGraph:
            if tet in includedTets:
                #all possible neighbors
                neighbors=fullMeshGraph[tet]
                #print(neighbors)
                #make sure chirality is preserved
                #neighbors=self.sortNeighbors(tet, neighbors)
                adjacency[tet]=[]
                #print(neighbors)
                for neighbor in neighbors:
                    if neighbor in includedTets:
                        adjacency[tet].append(neighbor)
                    else:
                        adjacency[tet].append(-1)
        return adjacency

    def reorder(self, adjacency):
        start=list(adjacency.keys())[0]
        neighbors, points=self.sortNeighbors(start, adjacency[start])
        adjacency[start]=neighbors
        self.face_lookup={frozenset([0,1,2]):0, frozenset([0,1,3]):1, frozenset([0,3,2]):2, frozenset([1,2,3]):3}
        self.faceToFace=[0,3,2,1]
        visited={start, -1}
        for i,neighbor in enumerate(neighbors):
            if neighbor not in visited:
                if i==0:
                    prev_point_types={points[0]:0, points[1]:1, points[2]:2}
                elif i==1:
                    prev_point_types={points[0]:0, points[1]:1, points[3]:3}
                elif i==2:
                    prev_point_types={points[0]:0, points[2]:2, points[3]:3}
                else:
                    prev_point_types={points[1]:1, points[2]:2, points[3]:3}
                visited.add(neighbor)
                self.sub_reorder(neighbor, start, prev_point_types, frozenset(prev_point_types.values()), self.faceToFace[i], visited, adjacency)
        
        
    def sub_reorder(self, current, previous, prev_point_types, prev_type_set, spot, visited, adjacency):
        #pdb.set_trace()
        neighbors=adjacency[current]
        adjacency[current]=[-1,-1,-1,-1]
        adjacency[current][spot]=previous
        previous_points=self.tetrahedrons[previous]
        #print("current",self.tetrahedrons[current])
        #print("previous", previous_points)
        for neighbor in neighbors:
            if neighbor!=previous and neighbor!=-1:
                neighbor_points=self.tetrahedrons[neighbor]
                #print(neighbor_points)
                #2 points shared with previous
                intersection=neighbor_points.intersection(previous_points)
                #find the 1 point shared with current that is not shared with previous
                point3=neighbor_points.intersection(self.tetrahedrons[current]).difference(intersection).pop()
                point_types12=set(prev_point_types[x] for x in intersection)
                point_type3={0,1,2,3}.difference(prev_type_set)
                point_types=frozenset(point_types12.union(point_type3))
                #figure out which side it belongs on based on the 3 point types of the points it shares with current
                side=self.face_lookup[point_types]
                new_point_types={x:prev_point_types[x] for x in intersection}
                new_point_types[point3]=point_type3.pop()
                adjacency[current][side]=neighbor
                if neighbor not in visited:
                    visited.add(neighbor)
                    self.sub_reorder(neighbor, current, new_point_types, point_types, self.faceToFace[side],visited,adjacency)
                

    def sortNeighbors(self, current, neighbors):
        '''
        Sorts the neighbors of a tetrahedron to ensure that they are in a specific chiral order
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        current -- the id of the current tetrahedron
        neighbors -- list of its neighbors
        
        Returns: a sorted version of its neighbors that ensures that all neighbor lists are in the same chiral order.
        
        Sorts the neighbors of a tetrahedron to ensure that they are in a specific chiral order. Does so by selecting a base
        triangle on the tet, then finding a normal vector to that tet using the cross product. Then determines if the 4th point
        in the tet is on the same side of the triangle as the normal vector. If not, rearranges the order of the points in the 
        base triangle (which flips the direction of the normal). Then uses this point ordering to order the neighbors.
        '''
        names=[x for x in self.tetrahedrons[current]]
        point1=self.meshCoordinates[names[0]]
        point2=self.meshCoordinates[names[1]]
        point3=self.meshCoordinates[names[2]]
        point4=self.meshCoordinates[names[3]]
        
        #we begin by ordering all the points on the tetrahedron.
        # get vector normal to triangle formed by points1-3
        cross=np.cross(np.array([point2[0]-point1[0], point2[1]-point1[1], point2[2]-point1[2]]),
                 np.array([point3[0]-point1[0], point3[1]-point1[1], point3[2]-point1[2]]))
        #determine whether the normal vector is pointing toward or away from point 4. Reorder points 2 and 3 if necessary
        factor=((point1[0]-point4[0])*cross[0]+(point1[1]-point4[1])*cross[1]+(point1[2]-point4[2])*cross[2])/(cross[0]**2+cross[1]**2+cross[2]**2)
        if factor<0:
            temp=names[1]
            names[1]=names[2]
            names[2]=temp
        #we then use those ordered points to order the neighbors.
        neighbors2=[-1,-1,-1,-1]
        for neighbor in neighbors:
            points=self.tetrahedrons[neighbor]
            if names[0] in points:
                if names[1] in points:
                    if names[2] in points:
                        neighbors2[0]=neighbor
                    else:
                        neighbors2[1]=neighbor
                else:
                    neighbors2[2]=neighbor
            else:
                neighbors2[3]=neighbor
        return neighbors2, names

    def fixIndices(self, adjacency,maxIndex):
        '''
        Fixes the indices of the tetrahedrons in the adjacency list to eliminate gaps.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        adjacency - adjacency list of the final mesh.
        maxIndex -- the highest current index (pre re-indexing)
        
        Returns: the reindexed adjacency list
        '''
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

    def createFullMeshGraph(self):
        '''
        Returns an adjacency list for the original tesselation. Reads it from file if it exists, otherwise generates it (and saves it to file)
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Returns: the adjacency list for the original tesselation.
        '''
        #checks whether the meshing is running out of the current module or out of the Unite module, changes path to saved adjacency accordingly
        pyFile=inspect.stack()[-1].filename
        path="Tet_base_meshes/"
        if "Meshing" not in pyFile:
            path="mesh/"+path
        files=os.listdir(path)
        targetFile="full_mesh_graph_"+str(self.meshSize)+".bin"
        #if a saved version of the adjacency exists
        if targetFile in files:
            adjacency=dict()
            file=open(path+targetFile, "rb")
            #reads the number of tetrahedrons
            lengthArray=array('i')
            lengthArray.fromfile(file, 1)
            values=array('i')
            #reads the adjacency info
            values.fromfile(file,lengthArray[0]*5)
            #pdb.set_trace()
            #each tet has at most 4 neighbors, -1 is used as filler for those with less.
            for i in range(0,len(values),5):
                adjacency[values[i]]=[]
                for j in range(i+1, i+5):
                    if values[j]!=-1:
                        adjacency[values[i]].append(values[j]) 
        #if no saved version exists
        else:
            #create mapping from points to tetrahedrons they're in.
            pointsToTets=dict()
            for i,tet in enumerate(self.tetrahedrons):
                for point in tet:
                    if point not in pointsToTets:
                        pointsToTets[point]={i}
                    else:
                        pointsToTets[point].add(i)
            adjacency=dict()
            for i in range(len(self.tetrahedrons)):
                adjacency[i]=[]
            #pdb.set_trace()
            #identify neighbors based on shared points.
            for i, tet in enumerate(self.tetrahedrons):
                pointCombos=itertools.combinations(tet, 3)
                for combo in pointCombos:
                    #tets are neighbors if they share 3 points.
                    neighbors=(pointsToTets[combo[0]].intersection(pointsToTets[combo[1]])).intersection(pointsToTets[combo[2]])
                    for neighbor in neighbors:
                        if neighbor>i:
                            adjacency[neighbor].append(i)
                            adjacency[i].append(neighbor)
            #write to file for time savings in the future.
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

    def getAllTetPoints(self):
        '''
        Get all the internal points for all the tetrahedrons in the full tesselation.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Get all the internal points for all the tetrahedrons in the full tesselation. Tries to read it in from file,
        which saves a considerable amount of time. If there is no saved file, generates the points and then saves to 
        file for future use.
        '''
        #scale the tesselation to cover the shape
        scale=self.factor*self.span/self.meshSize
        #the theoretical midpoint of the mesh (in practice this needs to be adjusted to preserve symmetry)
        midMesh=(self.meshSize-1)/2.0
        #determines how far each point in the original tesselation needs to be moved to cover the shape
        displacements=((self.dimensions[1]+self.dimensions[0])/2-scale*(midMesh-.5),
                       (self.dimensions[3]+self.dimensions[2])/2-scale*midMesh,
                       (self.dimensions[5]+self.dimensions[4])/2-scale*midMesh)
        #checks whether the meshing is running out of the current module or out of the Unite module, changes path to saved file accordingly
        pyFile=inspect.stack()[-1].filename
        path="Tet_base_meshes/"
        if "Meshing" not in pyFile:
            path="mesh/"+path
        files=os.listdir(path)
        targetFile="base_mesh_"+str(self.meshSize)+".bin"
        #if the file exists
        if targetFile in files:
            intValues=array('i')
            file=open(path+targetFile, "rb")
            #numCoordinates, numTets, interiorPointsPerTet
            intValues.fromfile(file, 3)
            #read in coordinates
            coordinateArray=array('d')
            coordinateArray.fromfile(file, 3*intValues[0])
            #read in tetrahedrons
            tetArray=array('i')
            tetArray.fromfile(file, 4*intValues[1])
            #read in interior points
            interiorPointArray=array('d')
            interiorPointArray.fromfile(file, intValues[1]*intValues[2]*3)
            file.close()
            self.meshCoordinates=[]
            #breaks up read in coordinate values into tuples of 3
            for i in range(0, len(coordinateArray), 3):
                self.meshCoordinates.append((coordinateArray[i]*scale+displacements[0],
                                    coordinateArray[i+1]*scale+displacements[1],
                                    coordinateArray[i+2]*scale+displacements[2]))
            self.tetrahedrons=[]
            #breaks up read in tet poitns into lists of 4
            for i in range(0, len(tetArray),4):
                self.tetrahedrons.append({tetArray[i],tetArray[i+1], tetArray[i+2], tetArray[i+3]})
            self.tetPoints=dict()
            key=0
            #breaks up read in interior points into groups of 3, then assigns them to their tet.
            for i in range(0, len(interiorPointArray),3*intValues[2]):
                self.tetPoints[key]=[]
                for j in range(i, i+3*intValues[2],3):
                    self.tetPoints[key].append([interiorPointArray[j]*scale+displacements[0],
                                       interiorPointArray[j+1]*scale+displacements[1],
                                       interiorPointArray[j+2]*scale+displacements[2]])
                key+=1
        else:
            #determine the number of interior points to use
            if self.meshSize<5:
                self.accuracy=3
            elif self.meshSize<16:
                self.accuracy=2
            elif self.meshSize<26:
                self.accuracy=1
            else:
                self.accuracy=0
            #generate original irish tesselation
            self.tetrahedrons, old_coordinates, mesh_dimensions=createIrish([0,self.meshSize,0,self.meshSize,0,self.meshSize], self.numCores)
            print(mesh_dimensions)
            self.absTol=1e-7
            self.tetPoints=dict()
            self.meshCoordinates=[]
            #rescale for current shape
            for coordinate in old_coordinates:
                self.meshCoordinates.append((coordinate[0]*scale+displacements[0], coordinate[1]*scale+displacements[1], coordinate[2]*scale+displacements[2]))
            #generate interior points for each tet.
            for i,tet in enumerate(self.tetrahedrons):
                self.tetPoints[i]=self.generateTetPoints([old_coordinates[x] for x in tet])
            #store tesselation as binary file to not have to regenerate each time. Use binary for lossless save of floats. save unscaled version
            file=open(path+targetFile, "wb")
            numInteriorPoints=[1,5,24,63]
            intValues=array('i',[len(self.meshCoordinates), len(self.tetrahedrons), numInteriorPoints[self.accuracy]])
            intValues.tofile(file)
            for coordinate in old_coordinates:
                cArray=array('d',coordinate)
                cArray.tofile(file)
            for tet in self.tetrahedrons:
                tArray=array('i', tet)
                tArray.tofile(file)
            for key in range(len(self.tetrahedrons)):
                for point in self.tetPoints[key]:
                    pArray=array('d',point)
                    pArray.tofile(file)
            file.close()
            for key in self.tetPoints:
                for point in self.tetPoints[key]:
                    for i in range(3):
                        point[i]=point[i]*scale+displacements[i]
            print(len(self.tetPoints),len(self.meshCoordinates))

    def generateTetPoints(self, vertices):
        '''
        Get all the internal points for a tetrahedron given its vertices.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        vertices -- The coordinates of the vertices of the triangle
        
        Get all the internal points for a tetrahedron given its vertices. Finds points relative to the centroid, and can do
        different numbers of points per tet.
        '''
        #4 of the midpoints of the edges of the tet
        mid1=[(vertices[0][0]+vertices[1][0])/2,(vertices[0][1]+vertices[1][1])/2,(vertices[0][2]+vertices[1][2])/2]
        mid2=[(vertices[2][0]+vertices[3][0])/2,(vertices[2][1]+vertices[3][1])/2,(vertices[2][2]+vertices[3][2])/2]
        mid3=[(vertices[0][0]+vertices[2][0])/2,(vertices[0][1]+vertices[2][1])/2,(vertices[0][2]+vertices[2][2])/2]
        mid4=[(vertices[1][0]+vertices[3][0])/2,(vertices[1][1]+vertices[3][1])/2,(vertices[1][2]+vertices[3][2])/2]
        #find 2 medians of the tetrahedron
        v1=[mid1[0]-mid2[0], mid1[1]-mid2[1], mid1[2]-mid2[2]]
        v2=[mid3[0]-mid4[0], mid3[1]-mid4[1], mid3[2]-mid4[2]]
        #solve for their intersection, which is the centroid
        mat=np.array([[v1[0],-v2[0]],[v1[1],-v2[1]]])
        if not math.isclose(np.linalg.det(mat),0.0,rel_tol=1e-06, abs_tol=self.absTol):
            solution=np.array([mid4[0]-mid2[0],mid4[1]-mid2[1]])
        else:
            mat=np.array([[v1[0],-v2[0]],[v1[2],-v2[2]]])
            if not math.isclose(np.linalg.det(mat),0.0,rel_tol=1e-06, abs_tol=self.absTol):
                solution=np.array([mid4[0]-mid2[0],mid4[2]-mid2[2]])
            else:
                mat=np.array([[v1[1],-v2[1]],[v1[2],-v2[2]]])
                solution=np.array([mid4[1]-mid2[1],mid4[2]-mid2[2]])
        values=np.linalg.solve(mat, solution)
        centroid=[mid2[0]+v1[0]*values[0],mid2[1]+v1[1]*values[0],mid2[2]+v1[2]*values[0]]
        #24 points
        if self.accuracy==2:
            points=[]
            mid5=[(vertices[0][0]+vertices[3][0])/2,(vertices[0][1]+vertices[3][1])/2,(vertices[0][2]+vertices[3][2])/2]
            mid6=[(vertices[2][0]+vertices[1][0])/2,(vertices[2][1]+vertices[1][1])/2,(vertices[2][2]+vertices[1][2])/2]
            #determine check points based on centroid. points between the centroid and the vertices and edge midpoints
            for exteriorPoint in vertices+[mid1,mid2,mid3,mid4, mid5, mid6]:
                points.append([(exteriorPoint[0]-centroid[0])*.5+centroid[0],(exteriorPoint[1]-centroid[1])*.5+centroid[1], (exteriorPoint[2]-centroid[2])*.5+centroid[2]])
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            #finds the centroids of the triangles of the tet
            triangle_centroids=[]
            triangle_centroids.append([(2/3)*mid1[0]+vertices[2][0]/3, (2/3)*mid1[1]+vertices[2][1]/3, (2/3)*mid1[2]+vertices[2][2]/3])
            triangle_centroids.append([(2/3)*mid2[0]+vertices[1][0]/3, (2/3)*mid2[1]+vertices[1][1]/3, (2/3)*mid2[2]+vertices[1][2]/3])
            triangle_centroids.append([(2/3)*mid4[0]+vertices[0][0]/3, (2/3)*mid4[1]+vertices[0][1]/3, (2/3)*mid4[2]+vertices[0][2]/3])
            triangle_centroids.append([(2/3)*mid3[0]+vertices[3][0]/3, (2/3)*mid3[1]+vertices[3][1]/3, (2/3)*mid3[2]+vertices[3][2]/3])
            #determine check points based on centroid and triangle centroids.
            for exteriorPoint in triangle_centroids:
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            return points
        #1 point
        elif self.accuracy==0:
            return [centroid]
        #5 points. Does corners and centroid.
        elif self.accuracy==1:
            points=[]
            for exteriorPoint in vertices:
                points.append([(exteriorPoint[0]-centroid[0])*.9+centroid[0],(exteriorPoint[1]-centroid[1])*.9+centroid[1], (exteriorPoint[2]-centroid[2])*.9+centroid[2]])
            points.append(centroid)
            return points
        #63 points. same points as in accuracy level 2, plus centroid and extra points between centroid and vertices, centroid and midpoints, centroid and 
        #triangle centroids, and centroids of the sub triangles on each triangle (if each triangle is divided into 4 a la fischer ski logo)
        elif self.accuracy==3:
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
    
    def makePly(self, ply_path, name, tetrahedrons, coordinates, includedTets, spread=False):
        '''
        Makes a ply file of the mesh for visualization.
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        ply_path -- the path to write the file to
        name -- name to give the shape
        tetrahedrons -- list of tetrahedrons, stored as lists of 4 point ids
        coordinates -- the points of the mesh
        includedTets -- the tetrahedrons that are included in the mesh.
        spread -- boolean of whether the tets in the ply file should be spread out (so you can see individual tets) or adjacent to each other.
        '''
        #find all the points and tets in the mesh.
        tets2=[]
        coordinates2=[]
        includedPoints=set()
        for i, tet in enumerate(tetrahedrons):
            if i in includedTets:
                tets2.append(tet)
                for point in tet:
                    includedPoints.add(point)
        mapping=dict()
        #reindex the points.
        for index in includedPoints:
            mapping[index]=len(coordinates2)
            coordinates2.append(coordinates[index])
        triangles=set()
        triangleLookup=[]
        #spreads out points from the approximate center of mass of the shape
        if spread:
            coordinates, tetrahedrons=self.spreadCoordinates(tetrahedrons, coordinates)
        #finds triangles from the tets
        for old_info in tets2:
            info=[mapping[x] for x in old_info]
            info.sort()
            t=(info[0], info[1], info[2])
            if t not in triangles:
                triangles.add(t)
                triangleLookup.append(t)
            t=(info[0], info[1], info[3])
            if t not in triangles:
                triangles.add(t)
                triangleLookup.append(t)
            t=(info[0], info[2], info[3])
            if t not in triangles:
                triangles.add(t)
                triangleLookup.append(t)
            t=(info[1], info[2], info[3])
            if t not in triangles:
                triangles.add(t)
                triangleLookup.append(t)
        triangles=triangleLookup
        #create the ply file using the triangles and the points.
        self.generatePly(coordinates2, triangles, ply_path, name, spread)
    
    def generatePly(self, coordinates, triangles, file_path, name, spread):
        '''
        Creates the actual plyfile 
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        file_path -- the path to write the file to
        name -- name to give the shape
        coordinates -- the coordinates for the file
        triangles -- the triangles for the file
        includedTets -- the tetrahedrons that are included in the mesh.
        spread -- boolean of whether the tets in the ply file should be spread out (so you can see individual tets) or adjacent to each other.
        '''
        pointList=[]
        faceList=[]
        #randomly color the tets if they are spread
        if spread:
            for i,p in enumerate(coordinates):
                if i%4==0:
                    color=[random.randint(0,255),random.randint(0,255),random.randint(0,255)]
                    colorSum=color[0]+color[1]+color[2]
                    if colorSum<200:
                        factor=200/colorSum
                        color=[math.ceil(color[0]*factor), math.ceil(color[1]*factor), math.ceil(color[2]*factor)]
                self.point(p[0],p[1],p[2], color, pointList)
        else:
            for p in coordinates:
                self.point(p[0],p[1],p[2], (250,150,150), pointList)
        for t in triangles:
            self.face3(t[0],t[1],t[2],faceList)
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
    
    def face3(self, a,b,c,facelist):
        '''
        Makes a 3 cornered face in the ply format using provided point ids
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        a -- the first point id
        b -- the second point id
        c -- the third point id
        facelist -- the list of face codes for the ply file
        '''
        facelist.append("3 "+str(a)+" "+str(b)+" "+str(c)+" \n")

    def spreadCoordinates(self, tetrahedrons, coordinates):
        '''
        Spreads out tetrahedrons from the approximate center of mass of the entire shape
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        tetrahedrons -- the tets in the final mesh
        coordinates -- the coordinates of the final mesh
        
        Returns: the spread out coordinates and tets referencing those coordinates.
        '''
        #find approximate center of mass of shape
        com=[0.0,0.0,0.0]
        total_vol=0.0
        radius=0.0
        com_ts=[]
        for i,tet in enumerate(tetrahedrons):
            point1=coordinates[tet[0]]
            point2=coordinates[tet[1]]
            point3=coordinates[tet[2]]
            point4=coordinates[tet[3]]
            com_t=([(point1[0]+point2[0]+point3[0]+point4[0])/4.0,(point1[1]+point2[1]+point3[1]+point4[1])/4.0, (point1[2]+point2[2]+point3[2]+point4[2])/4.0])
            com_ts.append(com_t)
            vol=self.volume(tet, coordinates)
            total_vol+=vol
            for i in range(3):
                com[i]+=com_t[i]*vol
        for i in range(3):
            com[i]=com[i]/total_vol
        #find distance from com
        for p in coordinates:
            d=self.dist(com,p)
            if d>radius:
                radius=d
        new_coordinates=[]
        new_tetrahedrons=[]
        for i,tet in enumerate(tetrahedrons):
            point1=coordinates[tet[0]]
            point2=coordinates[tet[1]]
            point3=coordinates[tet[2]]
            point4=coordinates[tet[3]]
            com_t=com_ts[i]
            d=self.dist(com, com_t)
            scale=3*d/radius
            #print(scale)
            #spread out points
            newPoint1=[point1[0]+scale*(com_t[0]-com[0]), point1[1]+scale*(com_t[1]-com[1]), point1[2]+scale*(com_t[2]-com[2])]
            newPoint2=[point2[0]+scale*(com_t[0]-com[0]), point2[1]+scale*(com_t[1]-com[1]), point2[2]+scale*(com_t[2]-com[2])]
            newPoint3=[point3[0]+scale*(com_t[0]-com[0]), point3[1]+scale*(com_t[1]-com[1]), point3[2]+scale*(com_t[2]-com[2])]
            newPoint4=[point4[0]+scale*(com_t[0]-com[0]), point4[1]+scale*(com_t[1]-com[1]), point4[2]+scale*(com_t[2]-com[2])]
            #print(point1, newPoint1)
            new_coordinates.append(newPoint1)
            new_coordinates.append(newPoint2)
            new_coordinates.append(newPoint3)
            new_coordinates.append(newPoint4)
            new_tetrahedrons.append([len(new_coordinates)-4,len(new_coordinates)-3,len(new_coordinates)-2,len(new_coordinates)-1])
        return new_coordinates, new_tetrahedrons
    
    def dist(self, point1, point2):
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

    def volume(self, points, coordinates):
        '''
        Calculates the volume of a tetrahedron
        
        Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 3/28/19
        
        Arguments:
        points -- the ids of the vertices of the tet
        coordinates -- the coordinates of the mesh
        
        Returns: the volume
        '''
        U=self.dist(coordinates[points[0]], coordinates[points[1]])
        V=self.dist(coordinates[points[1]], coordinates[points[2]])
        W=self.dist(coordinates[points[2]], coordinates[points[0]])
        u=self.dist(coordinates[points[2]], coordinates[points[3]])
        v=self.dist(coordinates[points[0]], coordinates[points[3]])
        w=self.dist(coordinates[points[1]], coordinates[points[3]])
        X=(w-U+v)*(U+v+w)
        x=(U-v+w)*(v-w+U)
        Y=(u-V+w)*(V+w+u)
        y=(V-w+u)*(w-u+V)
        Z=(v-W+u)*(W+u+v)
        z=(W-u+v)*(u-v+W)
        a=(x*Y*Z)**0.5
        b=(y*Z*X)**0.5
        c=(z*X*Y)**0.5
        d=(x*y*z)**0.5
        return (((-a+b+c+d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d))**0.5)/(192*u*v*w)

if __name__=="__main__":
    if 1==1:
        tm=TetMeshing()
        tm.run("examples/simpleframe.STL", "simpleframe", 256, 4)
        #print(tm.run("examples/sphere1.stl", "sphere1", 32, 4))
    if 0==1:
        hm=HexMeshing()
        hm.run("examples/simpleframe.STL", "simpleframe", 128, 4)
    if 0==1:
        tm=TetMeshing()
        tm.run("examples/lower_d.stl", "lower_d_tet", 1024,4)