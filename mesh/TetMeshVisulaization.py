import re, random, time, math
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.path as mpath
import matplotlib.patches as patches
#from IPython import get_ipython
import pdb
#get_ipython().run_line_magic('matplotlib', 'qt')
def dist(point1, point2):
    return ((point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2)**0.5

def volume(points, coordinates):
    U=dist(coordinates[points[0]], coordinates[points[1]])
    V=dist(coordinates[points[1]], coordinates[points[2]])
    W=dist(coordinates[points[2]], coordinates[points[0]])
    u=dist(coordinates[points[2]], coordinates[points[3]])
    v=dist(coordinates[points[0]], coordinates[points[3]])
    w=dist(coordinates[points[1]], coordinates[points[3]])
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

def spreadCoordinates(coordinates, tetrahedrons, volumes):
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
        vol=volumes[i]
        total_vol+=vol
        for i in range(3):
            com[i]+=com_t[i]*vol
    for i in range(3):
        com[i]=com[i]/total_vol
    for p in coordinates[1:]:
        d=dist(com,p)
        if d>radius:
            radius=d
    new_coordinates=[[]]
    new_tetrahedrons=[]
    for i,tet in enumerate(tetrahedrons):
        point1=coordinates[tet[0]]
        point2=coordinates[tet[1]]
        point3=coordinates[tet[2]]
        point4=coordinates[tet[3]]
        com_t=com_ts[i]
        d=dist(com, com_t)
        scale=3*d/radius
        #print(scale)
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
    

def plotCgal(meshFile, name, spread=False, hist=False):
    mesh=open(meshFile)
    triangles=set()
    triangleLookup=[]
    tetrahedrons=[]
    coordinates=[[]]
    minv=-1
    maxv=-1
    mint=[]
    maxt=[]
    numPoints=0
    volumes=[]
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
        if 0 in info:
            pdb.set_trace()
        vol=volume(info, coordinates)
        volumes.append(vol)
        if vol>maxv:
            maxv=vol
            maxt=info
        elif vol<minv or minv==-1:
            minv=vol
            mint=info
    if spread:
        coordinates, tetrahedrons=spreadCoordinates(coordinates, tetrahedrons, volumes)
    for info in tetrahedrons:
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
    makePly(coordinates, triangles, meshFile[:-5]+"_ply.ply", name, spread)
    if hist:
        fig, ax = plt.subplots()
        ax.hist(volumes, 100, density=1)
        plt.show()
    mesh.close()


def makePly(coordinates, triangles, file_path, name, spread):
    pointList=[]
    faceList=[]
    if spread:
        for i,p in enumerate(coordinates[1:]):
            if i%4==0:
                color=[random.randint(0,255),random.randint(0,255),random.randint(0,255)]
                colorSum=color[0]+color[1]+color[2]
                if colorSum<200:
                    factor=200/colorSum
                    color=[math.ceil(color[0]*factor), math.ceil(color[1]*factor), math.ceil(color[2]*factor)]
            point(p[0],p[1],p[2], color, pointList)
    else:
        for p in coordinates[1:]:
            point(p[0],p[1],p[2], (150,150,150), pointList)
    for t in triangles:
        face3(t[0],t[1],t[2],faceList)
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

def point(x,y,z,color,pointlist):
    pointlist.append(str(x)+" "+str(y)+" "+str(z)+" "+str(color[0])+" "+str(color[1])+" "+str(color[2])+" \n")

def face3(a,b,c,facelist):
    facelist.append("3 "+str(a-1)+" "+str(b-1)+" "+str(c-1)+" \n")

def plotTriangles(meshFile, gradient=False):
    mesh=open(meshFile)
    count=-1
    maxCount=-1
    reading1=False
    reading2=False
    triangles=set()
    triangleLookup=[]
    coordinates=[]
    minv=-1
    maxv=-1
    mint=[]
    maxt=[]
    #interpret Mesh file output from gmsh
    for line in mesh:
        #read triangles
        if reading2:
            count+=1
            if count==0:
                maxCount=int(line[:-1])
            elif count>0:
                spaces=[m.start() for m in re.finditer(" ", line)]
                if line[spaces[0]+1:spaces[0]+8]=="2 2 0 1":
                    info=[int(line[spaces[-3]+1:spaces[-2]]),int(line[spaces[-2]+1:spaces[-1]]),int(line[spaces[-1]+1:])]
                    info.sort()
                    info=tuple(info)
                    if info not in triangles:
                        triangles.add(info)
                        triangleLookup.append(info)
                elif line[spaces[0]+1:spaces[0]+8]=="4 2 0 1":
                    info=[int(line[spaces[-4]+1:spaces[-3]]), int(line[spaces[-3]+1:spaces[-2]]),int(line[spaces[-2]+1:spaces[-1]]),int(line[spaces[-1]+1:])]
                    vol=volume(info, coordinates)
                    if vol>maxv:
                        maxv=vol
                        maxt=info
                    elif vol<minv or minv==-1:
                        minv=vol
                        mint=info
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
    triangles=triangleLookup
    triangles=sortTriangles(triangles, coordinates)
    convertCoordinates(coordinates)
    figure=plt.figure(figsize=(10,10), frameon=False)
    axes=figure.add_axes([0,0,1,1])
    axes.set_axis_off()
    xmin=coordinates[0][0]
    ymin=coordinates[0][1]
    xmax=coordinates[0][0]
    ymax=coordinates[0][1]
    for coordinate in coordinates:
        if coordinate[0]<xmin:
            xmin=coordinate[0]
        elif coordinate[0]>xmax:
            xmax=coordinate[0]
        if coordinate[1]<ymin:
            ymin=coordinate[1]
        elif coordinate[1]>ymax:
            ymax=coordinate[1]
    axes.set_xbound(xmin,xmax)
    axes.set_ybound(ymin,ymax)
    plt.show()
    interval=60/len(triangles)
    count=0
    red=.2
    blue=.2
    green=.2
    colorIncrement=2.4/len(triangles)
    print(len(triangles))
    print(minv,maxv)
    print(mint,maxt)
    for triangle in triangles:
        if gradient:
            if red>=1:
                red=1
                if blue>=1:
                    blue=1
                    if green>=1:
                        green=1
                    else:
                        green+=colorIncrement
                        if green>1:
                            green=1
                else:
                    blue+=colorIncrement
                    if blue>1:
                        blue=1
            else:
                red+=colorIncrement
                if red>1:
                    red=1
        else:            
            red=random.uniform(.2,1)
            green=random.uniform(.2,1)
            blue=random.uniform(.2,1)
        path=[]
        path.append((Path.MOVETO,(coordinates[triangle[0]-1][0],coordinates[triangle[0]-1][1])))
        path.append((Path.LINETO,(coordinates[triangle[1]-1][0],coordinates[triangle[1]-1][1])))
        path.append((Path.LINETO,(coordinates[triangle[2]-1][0],coordinates[triangle[2]-1][1])))
        path.append((Path.LINETO,(coordinates[triangle[0]-1][0],coordinates[triangle[0]-1][1])))
        codes, verts = zip(*path)
        path = mpath.Path(verts, codes)
        patch = patches.PathPatch(path, linewidth=2, fill=True, aa=True, capstyle="round", facecolor=(red, green, blue), edgecolor=(red-.2, green-.2, blue-.2))
        axes.add_patch(patch)
        
        if len(triangles)<1000:
            time.sleep(interval)
            figure.canvas.draw()
            figure.canvas.flush_events()
        elif count%200==0:
            figure.canvas.draw()
            figure.canvas.flush_events()
        count+=1
        if count%500==0:
            print(count)
    figure.canvas.draw()
    figure.canvas.flush_events()
def convertCoordinates(coordinates):
    yadjust=.5
    xadjust=.5
    for coordinate in coordinates:
        coordinate[0]=coordinate[0]-xadjust*coordinate[2]
        coordinate[1]=coordinate[1]-yadjust*coordinate[2]

def sortTriangles(triangles, coordinates):
    trianglesc=[]
    for i in range(len(triangles)):
        triangle=triangles[i]
        d=coordinates[triangle[0]-1][0]+coordinates[triangle[0]-1][1]+coordinates[triangle[0]-1][2]
        for j in range(1,3):
            d2=coordinates[triangle[j]-1][0]+coordinates[triangle[j]-1][1]+coordinates[triangle[j]-1][2]
            if d2>d:
                d=d2
        trianglesc.append((d,i))
    trianglesc.sort()
    triangles2=[]
    for t in trianglesc:
        triangles2.append(triangles[t[1]])
    return triangles2
        
def sortTriangles2(triangles, coordinates):
    trianglesc=[]
    for i in range(len(triangles)):
        #trianglesc.append(((coordinates[triangles[i][0]-1][0]+coordinates[triangles[i][1]-1][0]+coordinates[triangles[i][2]-1][0])/3.0,
        #                   (coordinates[triangles[i][0]-1][1]+coordinates[triangles[i][1]-1][1]+coordinates[triangles[i][2]-1][1])/3.0,
        #                   (coordinates[triangles[i][0]-1][2]+coordinates[triangles[i][1]-1][2]+coordinates[triangles[i][2]-1][2])/3.0, i))
        trianglesc.append((max(coordinates[triangles[i][0]-1][0],coordinates[triangles[i][1]-1][0],coordinates[triangles[i][2]-1][0]),
                           max(coordinates[triangles[i][0]-1][1],coordinates[triangles[i][1]-1][1],coordinates[triangles[i][2]-1][1]),
                           min(coordinates[triangles[i][0]-1][2],coordinates[triangles[i][1]-1][2],coordinates[triangles[i][2]-1][2]), i))
    trianglesc=sort(trianglesc, 0, True)
    trianglesc=sort(trianglesc, 1, True)
    trianglesc=sort(trianglesc, 2, False)
    triangles2=[]
    for t in trianglesc:
        triangles2.append(triangles[t[3]])
    return triangles2
def sort(triangles, field, ascending):
    if len(triangles)>1:
        a=sort(triangles[:len(triangles)//2], field, ascending)
        b=sort(triangles[len(triangles)//2:], field, ascending)
        c=[]
        i=0
        j=0
        while i<len(a) or j<len(b):
            if i<len(a):
                if j<len(b):
                    if ascending:
                        if a[i][field]<=b[j][field]:
                            c.append(a[i])
                            i+=1
                        else:
                            c.append(b[j])
                            j+=1
                    else:
                        if a[i][field]>=b[j][field]:
                            c.append(a[i])
                            i+=1
                        else:
                            c.append(b[j])
                            j+=1
                else:
                    c=c+a[i:]
                    return c
            else:
                c=c+b[j:]
                return c
        return c
    else:
        return triangles
if __name__=="__main__":
    if 0==1:
        plotTriangles("examples/simpleframe_tetmesh0.msh")
    if 0==1:
        plotTriangles("examples/coolring2_tetmesh0.msh")
    if 0==1:
        plotTriangles("examples/coolring2_tetmesh0.msh", True)
    if 0==1:
        plotTriangles("examples/sphere.msh", True)
    if 0==1:
        plotCgal("examples/out.mesh", "cgal_test_sphere2", True)
    if 1==1:
        plotCgal("examples/test_cut_out.mesh", "test_cut_out", False)
    if 0==1:
        points=[0,1,2,3]
        coordinates=[[0,0,0],[-.25,.5,0],[.25,.5,0],[0,.25,.5]]
        print(volume(points, coordinates))
        