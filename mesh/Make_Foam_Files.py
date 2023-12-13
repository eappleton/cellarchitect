
# coding: utf-8

# In[1]:

import numpy as np
import os
import shutil
import pdb
import math
import subprocess
import binascii
import inspect
from OpenFoam import interpretFoam, downSize, fullyConnected


# In[2]:

def readSTL(stlPath):
    '''
    Helper function for makeFiles that reads an STL file and returns its extreme xyz coordinates and its triangles.
    
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
        numtri=interpretInt(stl2.read(4))
        #loops through all the triangles
        for i in range(numtri):
            stl2.read(12)
            triangle=[]
            for j in range(3):
                x=interpretFloat(stl2.read(4))
                y=interpretFloat(stl2.read(4))
                z=interpretFloat(stl2.read(4))
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


# In[3]:

def fixSTL(stlPath):
    '''
    Helper function for readStl. Fixes illegal headers in binary stl files.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/16/2017
    
    Keyword arguments:
    stlPath -- String that contains the path to a binary STL file

    Copies a binary STL file but rewrites its first 80 bits to 00. This is necessary because some programs that create
    binary STL files violate convention and place the word "solid" in the first 80 bits, which causes openfoam to
    believe that the file is ascii.
    '''
    stl=open(stlPath,"rb")
    path=stlPath
    while os.path.exists(path):
        path=path[:-4]+"1.stl"
    stl2=open(path, "wb")
    stl.read(80)
    for i in range(80):
        stl2.write(binascii.unhexlify("00"))
    temp=stl.read(2)
    while temp:
        stl2.write(temp)
        temp=stl.read(2)
    stl.close()
    stl2.close()
    return path


# In[4]:

def interpretFloat(fourbytes):
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


# In[5]:

def interpretInt(fourbytes):
    '''
    Helper function for readSTL that takes a 4 byte positive int read from an STL file and returns a python int.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 7/31/2017
    
    Keyword arguments:
    fourbytes -- bytearray of length 4 read from an STL file that represents a positive integer

    Takes a bytearray of length 4 read from an STL file and returns a python int. The bytes are in reverse order.
    '''
    return fourbytes[0]+(fourbytes[1]*256)+(fourbytes[2]*65536)+(fourbytes[3]*16777216)


# In[6]:

def blockMeshDict(dimensions, blocksize, filedest, templatepath):
    #pdb.set_trace()
    '''
    Helper function for makeFiles that creates a blockMeshDict for OpenFoam from a template.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/1/2017
    
    Keyword arguments:
    dimensions -- list containing the extreme xyz coordinates of the STL shape.
    blocksize -- the edge length of the cubes to be used in the mesh.
    filedest -- The path to write the resulting file to.
    templatepath -- The path of the template file.

    Creates a blockMeshDict for the OpenFoam blockMesh process based on the given dimensions and blocksize.
    '''
    template=open(templatepath)
    output=open(filedest,"w")
    for i in range(20):
        output.write(template.readline())
    #finds the center of the rectangular solid
    midx=(dimensions[0]+dimensions[1])/2
    midy=(dimensions[2]+dimensions[3])/2
    midz=(dimensions[4]+dimensions[5])/2
    print("mids",midx,midy,midz)
    print("blocksize",blocksize)
    #determines the number of blocks needed in each direction by calculing the number of blocks that will
    #fit between the center and the maximum for each coordinate, then adding 1 (thus the original shape is
    #entirely contained by the blockMesh), and multiplying by 2 (since we also need blocks between the
    #center and the minimum).
    numx=((dimensions[1]-midx)//blocksize+1)*2
    numy=((dimensions[3]-midy)//blocksize+1)*2
    numz=((dimensions[5]-midz)//blocksize+1)*2
    #updates the dimensions to be the dimensions of the blockMesh (as opposed to the dimensions of the
    #original shape)
    dimensions[0]=-(numx/2)*blocksize+midx
    dimensions[1]=(numx/2)*blocksize+midx
    dimensions[2]=-(numy/2)*blocksize+midy
    dimensions[3]=(numy/2)*blocksize+midy
    dimensions[4]=-(numz/2)*blocksize+midz
    dimensions[5]=(numz/2)*blocksize+midz
    line="    ( "+str(dimensions[0])+" "+str(dimensions[2])+" "+str(dimensions[4])+")\n"
    output.write(line)
    line="    ( "+str(dimensions[1])+" "+str(dimensions[2])+" "+str(dimensions[4])+")\n"
    output.write(line)
    line="    ( "+str(dimensions[1])+" "+str(dimensions[3])+" "+str(dimensions[4])+")\n"
    output.write(line)
    line="    ( "+str(dimensions[0])+" "+str(dimensions[3])+" "+str(dimensions[4])+")\n"
    output.write(line)
    line="    ( "+str(dimensions[0])+" "+str(dimensions[2])+" "+str(dimensions[5])+")\n"
    output.write(line)
    line="    ( "+str(dimensions[1])+" "+str(dimensions[2])+" "+str(dimensions[5])+")\n"
    output.write(line)
    line="    ( "+str(dimensions[1])+" "+str(dimensions[3])+" "+str(dimensions[5])+")\n"
    output.write(line)
    line="    ( "+str(dimensions[0])+" "+str(dimensions[3])+" "+str(dimensions[5])+")\n"
    output.write(line)
    for i in range(8):
        template.readline()
    for i in range(4):
        output.write(template.readline())
    line="    hex (0 1 2 3 4 5 6 7) ("+str(int(numx))+" "+str(int(numy))+" "+str(int(numz))+") simpleGrading (1 1 1)\n"
    output.write(line)
    template.readline()
    for i in range(25):
        output.write(template.readline())
    template.close()
    output.close()
    return dimensions


# In[7]:

def surfaceFeatureExtractDict(name, filedest, templatepath):
    '''
    Helper function for makeFiles that creates a surfaceFeatureExtractDict for OpenFoam from a template.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/1/2017
    
    Keyword arguments:
    name -- the name of an STL file
    filedest -- The path to write the resulting file to.
    templatepath -- The path of the template file.

    Creates a surfaceFeatureExtractDict for the OpenFoam surfaceFeatureExtract process.
    '''
    template=open(templatepath)
    output=open(filedest,"w")
    for i in range(16):
        output.write(template.readline())
    output.write(name+".stl\n")
    template.readline()
    for i in range(22):
        output.write(template.readline())
    template.close()
    output.close()


# In[8]:

def snappyHexMeshDict(name, filedest, templatepath, dimensions, triangles, blocksize):
    '''
    Helper function for makeFiles that creates a snappyHexMeshDict for OpenFoam from a template.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/1/2017
    
    Keyword arguments:
    name -- The name of an STL file
    filedest -- The path to write the resulting file to.
    templatepath -- The path of the template file.
    dimensions -- A list containing the dimensions of the blockMesh.
    triangles -- A list containing the triangles read from the STL file.
    blocksize -- The edge length of the blocks in the blockMesh.

    Creates a snappyHexMeshDict for the OpenFoam surfaceFeatureExtract process. Returns true if the
    dictionary creation was successful, false if there was an issue with finding an interior point.
    '''
    template=open(templatepath)
    output=open(filedest,"w")
    for i in range(31):
        output.write(template.readline())
    output.write("    "+name+".stl\n")
    template.readline()
    output.write(template.readline())
    output.write(template.readline())
    output.write("        name "+name+";\n")
    template.readline()
    for i in range(45):
        output.write(template.readline())
    output.write("            file \""+name+".eMesh\";\n")
    template.readline()
    for i in range(17):
        output.write(template.readline())
    output.write("        "+name+"\n")
    template.readline()
    for i in range(24):
        output.write(template.readline())
    output.write("        "+name+"\n")
    template.readline()
    for i in range(15):
        output.write(template.readline())
    #uses the findInteriorPoint function to find a point on the inside of the original STL shape.
    print("prepoint")
    point=findInteriorPoint(dimensions, triangles, blocksize)
    print("postpoint")
    #checks that an interior point was indeed found
    if point != None:
        output.write("    locationInMesh ("+str(point[0])+" "+str(point[1])+" "+str(point[2])+"); // Inside point\n")
        template.readline()
        template.readline()
        for i in range(171):
            output.write(template.readline())
        template.close()
        output.close()
        return True
    else:
        template.close()
        output.close()
        return False


# In[9]:

def findInteriorPoint(dimensions, triangles, blocksize, factor=1):
    '''
    Helper function for snappyHexMeshDict that finds a point interior to an STL shape.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/1/2017
    
    Keyword arguments:
    dimensions -- A list containing the dimensions of the blockMesh.
    triangles -- A list containing the triangles read from the STL file.
    blocksize -- The edge length of the blocks in the blockMesh.
    factor -- Number of extra times to divide the step size. Defaults to 1 which is a step size of 4xblocksize.

    Finds a point interior to an STL shape. Starts at the center of the blockMesh and checks points in cubic
    shells radiating outward.
    '''
    #pdb.set_trace()
    #determines the number of steps to take based on the factor and blocksize. With factor=1,
    #numx=the number of blocks in the x direction. This actually means that we are evaluating
    #steps that are 4x the blocksize, because numsteps is the number of steps between the
    #midpoint and the edge of the blockMesh.
    numx=(dimensions[1]-dimensions[0])*blocksize*factor/8
    numy=(dimensions[3]-dimensions[2])*blocksize*factor/8
    numz=(dimensions[5]-dimensions[4])*blocksize*factor/8
    numsteps=int(max(numx,numy,numz))
    #determines the center of the blockMesh
    midx=(dimensions[1]+dimensions[0])/2
    midy=(dimensions[3]+dimensions[2])/2
    midz=(dimensions[5]+dimensions[4])/2
    #loops through the steps away from the center
    for step in range(numsteps):
        #functionality for higher factors that ensures that points tested at previous lower factors will
        #not be repeated. Higher factors will only be called recursively and will always be half the
        #previous factor. Thus when factor>1, we can ignore all even steps because they were alread tested
        #at lower factors.
        if factor==1 or (step%2)==1:
            #Tests  points on the positive x face (constant x where x>=midx, variable y and z)
            #calculates x
            x=midx+(step*(blocksize*4/factor))
            #if x is not ouside of the blockMesh
            if x < dimensions[1]:
                #for y in [-step,step]
                for i in range(-step,step+1):
                    #calculates y
                    y=midy+(i*(blocksize*4/factor))
                    #if y is not ouside of the blockMesh
                    if y > dimensions[2] and y < dimensions[3]:
                        #for z in [-step,step]
                        for j in range(-step,step+1):
                            #calculates z
                            z=midz+(j*(blocksize*4/factor))
                            #if z is not ouside of the blockMesh
                            if z > dimensions[4] and z < dimensions[5]:
                                #checks the point
                                if (not onCell(x, y, z, dimensions, blocksize)) and checkPoint(x,y,z,triangles, blocksize):
                                    #if the point was inside, returns the point.
                                    return [x,y,z]
            #Tests  points on the negative x face (constant x where x<=midx, variable y and z)
            x=midx-(step*(blocksize*4/factor))
            if x > dimensions[0]:
                for i in range(-step,step+1):
                    y=midy+(i*(blocksize*4/factor))
                    if y > dimensions[2] and y < dimensions[3]:
                        for j in range(-step,step+1):
                            z=midz+(j*(blocksize*4/factor))
                            if z > dimensions[4] and z < dimensions[5]:
                                if (not onCell(x, y, z, dimensions, blocksize)) and checkPoint(x,y,z,triangles, blocksize):
                                    return [x,y,z]
            #Tests  points on the postive y face (constant y where y>=midy, variable x and z)
            y=midy+(step*(blocksize*4/factor))
            if y < dimensions[3]:
                #for x in [-step+1,step-1]. We decrease the range of x because we have already checked the
                #points where x=-step or x=step in the previous 2 loops.
                for i in range(-step+1,step):
                    x=midx+(i*(blocksize*4/factor))
                    if x > dimensions[0] and x < dimensions[1]:
                        for j in range(-step,step+1):
                            z=midz+(j*(blocksize*4/factor))
                            if z > dimensions[4] and z < dimensions[5]:
                                if (not onCell(x, y, z, dimensions, blocksize)) and checkPoint(x,y,z,triangles, blocksize):
                                    return [x,y,z]
            #Tests  points on the negative y face (constant y where y<=midy, variable x and z)
            y=midy-(step*(blocksize*4/factor))
            if y > dimensions[2]:
                #for x in [-step+1,step-1]. We decrease the range of x because we have already checked the
                #points where x=-step or x=step in the first 2 loops.
                for i in range(-step+1,step):
                    x=midx+(i*(blocksize*4/factor))
                    if x > dimensions[0] and x < dimensions[1]:
                        for j in range(-step,step+1):
                            z=midz+(j*(blocksize*4/factor))
                            if z > dimensions[4] and z < dimensions[5]:
                                if (not onCell(x, y, z, dimensions, blocksize)) and checkPoint(x,y,z,triangles, blocksize):
                                    return [x,y,z]
            #Tests  points on the positive z face (constant z where z>=midz, variable x and y)
            z=midz+(step*(blocksize*4/factor))
            if z < dimensions[5]:
                #for x in [-step+1,step-1]. We decrease the range of x because we have already checked the
                #points where x=-step or x=step in the first 2 loops.
                for i in range(-step+1,step):
                    x=midx+(i*(blocksize*4/factor))
                    if x > dimensions[0] and x < dimensions[1]:
                        #for y in [-step+1,step-1]. We decrease the range of y because we have already checked the
                        #points where y=-step or y=step in the 3rd and 4th loops loops.
                        for j in range(-step+1,step):
                            y=midy+(j*(blocksize*4/factor))
                            if y > dimensions[2] and y < dimensions[3]:
                                if (not onCell(x, y, z, dimensions, blocksize)) and checkPoint(x,y,z,triangles, blocksize):
                                    return [x,y,z]
            #Tests  points on the negative z face (constant z where z<=midz, variable x and y)
            z=midz-(step*(blocksize*4/factor))
            if z > dimensions[4]:
                #for x in [-step+1,step-1]. We decrease the range of x because we have already checked the
                #points where x=-step or x=step in the first 2 loops.
                for i in range(-step+1,step):
                    x=midx+(i*(blocksize*4/factor))
                    if x > dimensions[0] and x < dimensions[1]:
                        #for y in [-step+1,step-1]. We decrease the range of y because we have already checked the
                        #points where y=-step or y=step in the 3rd and 4th loops loops.
                        for j in range(-step+1,step):
                            y=midy+(j*(blocksize*4/factor))
                            if y > dimensions[2] and y < dimensions[3]:
                                if (not onCell(x, y, z, dimensions, blocksize)) and checkPoint(x,y,z,triangles, blocksize):
                                    return [x,y,z]
    #if the factor is small enough and we haven't yet found an interior point, we increase the factor
    #in order to decrease the step size and recursively call this function
    if factor < 64:
        return findInteriorPoint(dimensions, triangles, blocksize, factor*2)
    #if the factor is relatively large, we give up and assume that there is an issue with the STL
    #shape (such as it not being closed) and return None.
    return None


# In[ ]:

def onCell(x,y,z, dimensions, blocksize):
    tol=blocksize*(1e-09)
    temp=(x-dimensions[0])/blocksize
    if math.isclose(round(temp), temp, rel_tol=1e-06, abs_tol=tol):
        return True
    temp=(y-dimensions[2])/blocksize
    if math.isclose(round(temp), temp, rel_tol=1e-06, abs_tol=tol):
        return True
    temp=(z-dimensions[4])/blocksize
    if math.isclose(round(temp), temp, rel_tol=1e-06, abs_tol=tol):
        return True
    return False


# In[10]:

def checkPoint(x,y,z,triangles,blocksize):
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
    tol=blocksize*(1e-09)
    np.seterr(all='raise')
    count=0
    problems=[]
    for triangle in triangles:
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
            xclose=math.isclose(triangle[0][0], x, rel_tol=1e-06, abs_tol=tol) or math.isclose(triangle[1][0], x, rel_tol=1e-06, abs_tol=tol) or math.isclose(triangle[2][0], x, rel_tol=1e-06, abs_tol=tol)
            #test if the triangle has at least one vector that is in the positive x direction from
            #the point. If not we can safely ignore the triangle
            if triangle[0][0]>x or triangle[1][0]>x or triangle[2][0]>x or xclose:
                #test if the positive x vector is at least in the right range to cross the triangle
                if (triangle[0][1]>=y or triangle[1][1]>=y or triangle[2][1]>=y) and (triangle[0][1]<=y or triangle[1][1]<=y or triangle[2][1]<=y):
                    if (triangle[0][2]>=z or triangle[1][2]>=z or triangle[2][2]>=z) and (triangle[0][2]<=z or triangle[1][2]<=z or triangle[2][2]<=z):
                        #pdb.set_trace()
                        #finds a vector that is normal to the plane that the triangle ABC is in using the cross product of
                        #vector A-B and vector A-C
                        line1=[triangle[0][0]-triangle[1][0],triangle[0][1]-triangle[1][1],triangle[0][2]-triangle[1][2]]
                        line2=[triangle[0][0]-triangle[2][0],triangle[0][1]-triangle[2][1],triangle[0][2]-triangle[2][2]]
                        cross=np.cross(line1,line2)
                        solution=[]
                        solution.append((cross[0]*triangle[0][0])+(cross[1]*triangle[0][1])+(cross[2]*triangle[0][2]))
                        solution.append(y)
                        solution.append(z)
                        #test if the positive x vector from the point intersects the plane that the triangle is in
                        try:
                            point=np.linalg.solve(np.array([cross,[0,1,0],[0,0,1]]), np.array(solution))
                            done=False
                        except np.linalg.LinAlgError:
                            done=True
                        #if the vector intersects the plane, we see if the point of intersection is in the triangle.
                        #we then check to see if the point is in the positive x direction
                        if not done:
                            xclose=math.isclose(point[0], x, rel_tol=1e-06, abs_tol=tol)
                            if point[0]>x or xclose:
                                
                                #up to here
                                
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
                                    solution=[point[0]-triangle[0][0],point[1]-triangle[0][1]]
                                    intest=np.linalg.solve(np.array([line1,line2]),np.array(solution))
                                elif not math.isclose((line1[0]*line3[1]), (line1[1]*line3[0]), rel_tol=1e-06, abs_tol=tol):
                                    solution=[point[0]-triangle[0][0],point[2]-triangle[0][2]]
                                    intest=np.linalg.solve(np.array([line1,line3]),np.array(solution))
                                else:
                                    solution=[point[1]-triangle[0][1],point[2]-triangle[0][2]]
                                    intest=np.linalg.solve(np.array([line2,line3]),np.array(solution))
                                #if the point is in the triangle
                                if (intest[0]>0 or math.isclose(intest[0],0, rel_tol=1e-06, abs_tol=tol)) and (intest[1]>0 or math.isclose(intest[1],0, rel_tol=1e-06, abs_tol=tol)) and ((intest[0]+intest[1])<1 or math.isclose(intest[0]+intest[1],1, rel_tol=1e-06, abs_tol=tol)):
                                    #if it's in the triangle but the x distance between the interesection point and
                                    #the original point was too close, then we find a new point entirely
                                    if xclose:
                                        return False
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
                                    count+=1
    #if the count is odd then the point is inside the shape.
    if int(count)%2 == 1:
        return True
    return False


# In[11]:

def checkPath(path):
    '''
    Helper function for makeFiles that checks if a path is suitable for OpenFoam.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/1/2017
    
    Keyword arguments:
    path -- String containing the path in question.

    Checks if a path is suitable for OpenFoam by checking for illegal characters or illegal character positions.
    '''
    for x in path:
        #checks that the path contains only permitted characters (letters, numbers, or the 5 symbols)
        if not (x.isalpha() or x.isdecimal() or x in "_/.()"):
            return False
    #checks that no name within the path begins with a number.
    slash=path.find("/")
    while (slash+1)<len(path):
        x=path[slash+1]
        if not (x.isalpha() or x=="_"):
            return False
        slash=path.find("/",slash+1)
    return True


# In[12]:

def makeFiles(stl,templates, maxCells, blocksize=-1):
    '''
    Helper function for master that creates necessary files for OpenFoam based on an stl shape.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/1/2017
    
    Keyword arguments:
    stl -- String that contains the path to an STL file
    templates -- String that contains the path to the folder that contains the OpenFoam templates.
    maxCells -- The maximum number of cells to use

    Creates necessary files to run OpenFoam blockMesh, surfaceFeatureExtract, and snappyHexMesh routines.
    '''
    path=stl[:-4]+"/"
    origpath=stl[:-4]+"/"
    #finds a path that will not cause problems for OpenFoam
    while not checkPath(path):
        path=path[:-path[::-1].find("/",1)]
    #creates a temporary directory that we will run OpenFoam on
    path=path+"temp/"
    name="temp"
    #makes sure that the directory we are creating does not already exist.
    while os.path.exists(path):
        path=path[:-1]+"1/"
        name=name+"1"
    os.makedirs(path+"system/",0o777,True)
    os.makedirs(path+"constant/triSurface/",0o777,True)
    os.makedirs(path+"constant/polyMesh/",0o777,True)
    os.makedirs(path+"0/",0o777,True)
    triangles=readSTL(stl)
    print("read stl")
    tempstl=path+"constant/triSurface/"+name+".stl"
    if triangles[7]:
        shutil.move(triangles[7],tempstl)
    else:
        shutil.copyfile(stl,tempstl)
    dimensions=triangles[:6]
    print(dimensions)
    triangles=triangles[6]
    if blocksize<0:
        #we calculate blocksize based on target number of cells. Divide by 2 to implement strategy of using smaller blocks, then converting to larger to reduce gaps
        blocksize=((((dimensions[1]-dimensions[0])*(dimensions[3]-dimensions[2])*(dimensions[5]-dimensions[4]))/(maxCells/4.0))**(1.0/3.0))/2
        blocksize=round(blocksize, -math.floor(math.log10(blocksize))+2)
    #will change dimensions so that dimensions becomes the dimensions of the blockmesh
    gridDimensions=blockMeshDict(dimensions, blocksize, path+"constant/polyMesh/blockMeshDict",templates+"blockMeshDict")
    print("blockmeshdict")
    shutil.copyfile(templates+"controlDict",path+"system/controlDict")
    shutil.copyfile(templates+"decomposeParDict",path+"system/decomposeParDict")
    shutil.copyfile(templates+"fvSchemes",path+"system/fvSchemes")
    shutil.copyfile(templates+"fvSolution",path+"system/fvSolution")
    shutil.copyfile(templates+"meshQualityDict",path+"system/meshQualityDict")
    os.makedirs(path+"system/caseDicts/",0o777,True)
    shutil.copyfile(templates+"meshQualityDict",path+"system/caseDicts/meshQualityDict")
    triangles=triangles[:-1]
    surfaceFeatureExtractDict(name, path+"system/surfaceFeatureExtractDict", templates+"surfaceFeatureExtractDict")
    print("surfaceextract")
    fine=snappyHexMeshDict(name, path+"system/snappyHexMeshDict", templates+"snappyHexMeshDict", dimensions, triangles, blocksize)
    print("snappyhex")
    if fine:
        return path, origpath, blocksize, gridDimensions
    else:
        print("Invalid STL shape")
        return None,None, None, None


# In[13]:

def runFoam(path):
    '''
    Helper function for master that sends commands to a linux shell that run openfoam
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/16/2017
    
    Keyword arguments:
    path -- The string path to the directory we wish to run openfoam in.
    '''
    os.system("cd "+path+"&& blockMesh")
    os.system("cd "+path+"&& surfaceFeatureExtract")
    os.system("cd "+path+"&& snappyHexMesh")


# In[14]:

def cleanup(path, origpath):
    '''
    Helper function for master that deletes files created while running openfoam.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/16/2017
    
    Keyword arguments:
    path - The path to the directory that we ran openfoam out of.
    origpath - the path to the directory we will send the results to.

    Moves files from the directory in which we ran openfoam to the directory where we want our results.
    Adds a file to faciliate opening in Paraview.
    '''
    old=None
    if os.path.exists(origpath+"system/"):
        if not old:
            old=origpath+"Old/"
            while os.path.exists(old):
                old=old[:-1]+"1/"
            os.makedirs(old,0o777,True)
        shutil.move(origpath+"system/",old+"system/")
    if os.path.exists(origpath+"constant/"):
        if not old:
            old=origpath+"Old/"
            while os.path.exists(old):
                old=old[:-1]+"1/"
            os.makedirs(old,0o777,True)
        shutil.move(origpath+"constant/",old+"constant/")
    if os.path.exists(origpath+"0/"):
        if not old:
            old=origpath+"Old/"
            while os.path.exists(old):
                old=old[:-1]+"1/"
            os.makedirs(old,0o777,True)
        shutil.move(origpath+"0/",old+"0/")
    if os.path.exists(origpath+"1/"):
        if not old:
            old=origpath+"Old/"
            while os.path.exists(old):
                old=old[:-1]+"1/"
            os.makedirs(old,0o777,True)
        shutil.move(origpath+"1/",old+"1/")
    shutil.move(path+"system/",origpath+"system/")
    shutil.move(path+"constant/",origpath+"constant/")
    shutil.move(path+"0/",origpath+"0/")
    shutil.move(path+"1/",origpath+"1/")
    tempfile=open(origpath+origpath[-origpath[-2::-1].find("/")-1:-1]+".foam","w")
    tempfile.close()
    os.rmdir(path)


# In[15]:

def master(stl, maxCells):
    '''
    Automates the snappyHexMesh process.
    
    Author: Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
    
    Date: 8/16/2017
    
    Keyword arguments:
    stl -- String that contains the path to an STL file
    maxCells -- The maximum number of cells 

    Runs openfoam scripts blockMesh, surfaceFeatureExtract, and snappyHexMesh on a STL file, creating required
    dictionaries.
    '''
    pyFile=inspect.stack()[-1].filename
    if pyFile=="Make_Foam_Files.py":
        templates="Foam_Templates/"
    else:
        templates="mesh/Foam_Templates/"
    while True:
        path0,path1,blocksize, gridDimensions=makeFiles(stl, templates, maxCells)
        if path0:
            runFoam(path0)
            cleanup(path0,path1)
            path2=path1+"1/polyMesh/"
            adjacency, coordinates=interpretFoam(path2+"owner",path2+"neighbour",path2+"faces",path2+"points", "", blocksize)
            print(adjacency)
            adjacency=downSize(adjacency, coordinates, gridDimensions, blocksize)[0]
            print(adjacency)
            iterations=0
            while not fullyConnected(adjacency):
                if iterations==7:
                    raise Exception("Non Connected Shape")
                iterations+=1
                path0,path1,blocksize, gridDimensions=makeFiles(stl, templates, maxCells, blocksize/(2**iterations))
                runFoam(path0)
                cleanup(path0,path1)
                path2=path1+"1/polyMesh/"
                adjacency, coordinates=interpretFoam(path2+"owner",path2+"neighbour",path2+"faces",path2+"points", "", blocksize)
                for i in range(iterations+1):
                    adjacency, coordinates, gridDimensions=downSize(adjacency, coordinates, gridDimensions, blocksize*(2**i))
            numBlocks=len(adjacency)
            if numBlocks*8>maxCells:
                path0, path1, blocksize, gridDimensions=makeFiles(stl, templates, maxCells, blocksize*1.5)
            elif numBlocks*16<maxCells:
                path0, path1, blocksize, gridDimensions=makeFiles(stl, templates, maxCells, blocksize*.75)
            else:
                return path1, adjacency
        else:
            raise Exception("Hex Meshing Issue")


# In[17]:
if __name__=="__main__":
    if 0==1:
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/logo.stl",1)
    if 0==1:
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/test_cylinder_2.stl",15)
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/test_cylinder_1.stl",15)
    if 0==1:
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/coolring.STL",5)
    if 0==1:
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/coolring2.STL",25)
    if 0==1:
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/simpleframe.STL",3)
    if 0==1:
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubesphere.STL",20)
    if 0==1:
        master("/mnt/c/Users/Tristan/Documents/ELM_Project/makefoamfilestest/tubespherebig.STL",128)
    if 0==1:
        print(master("/mnt/c/Users/Tristan/Documents/CAD_bio/src/python/mesh/examples/coolring2.STL",128))
    if 0==1:
        print(master("/mnt/c/Users/Tristan/Documents/CAD_bio/src/python/mesh/examples/sphere1.stl",64))
    if 1==1:
        print(master("/mnt/c/Users/Tristan/Documents/CAD_bio/src/python/mesh/examples/simpleframe.STL",128))

