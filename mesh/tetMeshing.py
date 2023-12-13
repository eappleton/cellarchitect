from Make_Foam_Files import readSTL
from hexMeshing import HexMeshing
import os, re
import pdb
import math, random
import numpy as np

def stlToOFF(stl_path):
    data=readSTL(stl_path)
    triangles=data[6]
    points=dict()
    triangles2=[]
    pointLookup=[]
    for triangle in triangles:
        newTriangle=[]
        for point in triangle:
            point2=tuple(point)
            if point2 in points:
                key=points[point2]
            else:
                key=len(points)
                points[point2]=key
                pointLookup.append(point2)
            newTriangle.append(key)
        triangles2.append(newTriangle)
    off_path=stl_path[:-3]+"off"
    file=open(off_path,"w")
    file.write("OFF\n\n")
    file.write(str(len(points))+" "+str(len(triangles))+"\n")
    for point in pointLookup:
        file.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")
    for triangle in triangles2:
        file.write("3 "+str(triangle[0])+" "+str(triangle[1])+" "+str(triangle[2])+"\n")
    return data[:6], triangles

def makeCgalFile(off_path, dimensions):
    template=open("cgal_template.cpp")
    outFile=open("cgal_file.cpp","w")
    count=0
    span=max(dimensions[1]-dimensions[0], dimensions[3]-dimensions[2], dimensions[5]-dimensions[4])
    facet_size=str(span*2)
    cell_size=str(span/1)
    #facet_size=str(span/16)
    #cell_size=str(span/32)
    for line in template:
        count+=1
        if count==30:
            loc1=line.find("=")
            outFile.write(line[:loc1+2]+"\""+off_path+"\";\n")
        elif count==48:
            loc1=line.find("facet_size")
            outFile.write(line[:loc1+11]+facet_size+", facet_distance="+cell_size+",\n")
        elif count==49:
            loc1=line.find("cell_size")
            outFile.write(line[:loc1+10]+cell_size+");\n")
        else:
            outFile.write(line)
    template.close()
    outFile.close()

def checkPointRelaxed(x,y,z,triangles,blocksize, triangleGrid, tgFactor):
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

def trimMesh(mesh_path, stl_triangles, dimensions):
    mesh=open(mesh_path)
    tetrahedrons=[]
    coordinates=[[]]
    pointLines=[""]
    pointSet=set()
    tetSet=set()
    triangles=[]
    triangleSet=set()
    numPoints=0
    while True:
        line=mesh.readline()
        if "Vertices" in line:
            line=mesh.readline()
            numPoints=int(line)
            break
    for i in range(numPoints):
        line=mesh.readline()
        pointLines.append(line)
        spaces=[m.start() for m in re.finditer(" ", line)]
        coordinates.append([float(line[:spaces[0]]), float(line[spaces[0]+1:spaces[1]]), float(line[spaces[1]+1:spaces[2]])])
    line=mesh.readline()
    if "Triangles" not in line:
        print("Problem")
    numTriangles=int(mesh.readline())
    for i in range(numTriangles):
        line=mesh.readline()
        spaces=[m.start() for m in re.finditer(" ", line)]
        info=[int(line[:spaces[0]]), int(line[spaces[0]+1:spaces[1]]), int(line[spaces[1]+1:spaces[2]])]
        triangles.append(info)
    line=mesh.readline()
    if "Tetrahedra" not in line:
        print("Problem")
    numTets=int(mesh.readline())
    for i in range(numTets):
        line=mesh.readline()
        spaces=[m.start() for m in re.finditer(" ", line)]
        info=[int(line[:spaces[0]]), int(line[spaces[0]+1:spaces[1]]), int(line[spaces[1]+1:spaces[2]]), int(line[spaces[2]+1:spaces[3]])]
        tetrahedrons.append(info)
    inCoord=set()
    hm=HexMeshing()
    triangleGrid,tgFactor=hm.triageTriangles(stl_triangles, dimensions)
    fake_blocksize=max(dimensions[1]-dimensions[0], dimensions[3]-dimensions[2], dimensions[5]-dimensions[4])/10
    for i,point in enumerate(coordinates[1:]):
        #if 50<point[0]<70 and 25<point[1]<40 and point[2]<8:
            #pdb.set_trace()
        if checkPointRelaxed(point[0],point[1],point[2], stl_triangles, fake_blocksize, triangleGrid, tgFactor):
            inCoord.add(i+1)
        #else:
        #    print(point)
    for i,tet in enumerate(tetrahedrons):
        count=0
        for x in tet:
            if x in inCoord:
                count+=1
                break
        if count>0:
            tetSet.add(i)
            for x in tet:
                pointSet.add(x)
    for i, triangle in enumerate(triangles):
        allIn=True
        for x in triangle:
            if x not in pointSet:
                allIn=False
                break
        if allIn:
            triangleSet.add(i)
    count=0
    pointMapping=dict()
    for i in range(1,len(pointLines)):
        if i in pointSet:
            pointMapping[i]=i-count
        else:
            count+=1
    outfile=open(mesh_path[:-5]+"_edit.mesh","w")
    outfile.write("MeshVersionFormatted 1\nDimension 3\nVertices\n")
    outfile.write(str(len(pointSet))+"\n")
    for i,line in enumerate(pointLines):
        if i in pointSet:
            outfile.write(line)
    outfile.write("Triangles\n"+str(len(triangleSet))+"\n")
    for i, triangle in enumerate(triangles):
        if i in triangleSet:
            outfile.write(str(pointMapping[triangle[0]])+" "+str(pointMapping[triangle[1]])+" "+str(pointMapping[triangle[2]])+" 1\n")
    outfile.write("Tetrahedra\n"+str(len(tetSet))+"\n")
    for i, tet in enumerate(tetrahedrons):
        if i in tetSet:
            outfile.write(str(pointMapping[tet[0]])+" "+str(pointMapping[tet[1]])+" "+str(pointMapping[tet[2]])+" "+str(pointMapping[tet[3]])+" 1\n")
    outfile.write("End\n")
    outfile.close()

def runCgal(path):
    os.system("cd "+path+" && g++ cgal_file.cpp -lCGAL -lgmp -lmpfr")
    os.system("cd "+path+" && ./a.out")
    os.remove("cgal_file.cpp")
    
if __name__=="__main__":
    if 0==1:
        dimensions, stl_triangles=stlToOFF("examples/sphere1.stl")
        makeCgalFile("examples/sphere1.off", dimensions)
        runCgal("/mnt/c/Users/Tristan/Documents/CAD_bio/src/python/mesh/")
    if 0==1:
        dimensions, stl_triangles=stlToOFF("examples/simpleframe.STL")
        makeCgalFile("examples/simpleframe.off", dimensions)
        runCgal("/mnt/c/Users/Tristan/Documents/CAD_bio/src/python/mesh/")
        trimMesh("examples/out.mesh", stl_triangles, dimensions)
    if 1==1:
        dimensions, stl_triangles=stlToOFF("examples/cube1.stl")
        makeCgalFile("examples/cube1.off", dimensions)
        runCgal("/mnt/c/Users/Tristan/Documents/CAD_bio/src/python/mesh/")
                