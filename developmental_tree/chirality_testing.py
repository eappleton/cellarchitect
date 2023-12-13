import itertools, pdb

def tetChirality():
    points=[(-1,0,-.707), (1,0,-.707), (0,-1,.707), (0,1,.707)]
    permutations=itertools.permutations(points[1:])
    colors=[(0,0,0),(200,0,0),(0,200,0),(0,0,200)]
    file_path="tet_chirality"
    for j,perm in enumerate(permutations):
        pointList=[]
        faceList=[]
        current_points=[points[0]]+list(perm)
        print(current_points)
        for i,p in enumerate(current_points):
            point(p[0],p[1],p[2], colors[i], pointList)
        for i in range(4):
            face_points=[0,1,2,3]
            face_points.remove(i)
            face3(face_points[0], face_points[1], face_points[2], faceList)
        outfile=open(file_path+str(j)+".ply", "w")
        outfile.write("ply\n")
        outfile.write("format ascii 1.0\n")
        outfile.write("obj_info tet\n")
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

def rotate(colors):
    rotations=[]
    rotations.append(tuple(colors))
    for i in range(23):
        if i%4==3:
            colors=turn(colors)
        else:
            colors=roll(colors)
        rotations.append(colors)
    return rotations

def intelligentPermutations(colors):
    perm1=itertools.permutations(range(6))
    final_perm=[]
    group1={1,2,4}
    group2={3,5,6}
    group1_indices=[0,1,3]
    group2_indices=[2,4,5]
    for p in perm1:
        invalid=False
        for i in group1_indices:
            x=p[i]
            if x in group2:
                invalid=True
                break
        if not invalid:
            for i in group2_indices:
                x=p[i]
                if x in group1:
                    invalid=True
                    break
        if not invalid:
            these_colors=[colors[0]]
            for x in p:
                these_colors.append(colors[x+1])
            these_colors.append(colors[7])
            final_perm.append(these_colors)
    return final_perm

def roll(colors):
    return (colors[1],colors[5],colors[3],colors[7],colors[0], colors[4], colors[2], colors[6])

def turn(colors):
    return (colors[2], colors[3], colors[6], colors[7], colors[0], colors[1], colors[4], colors[5])

def hexChirality():
    points=[(-1,-1,-1),(-1,1,-1),(-1,-1,1),(-1,1,1),(1,-1,-1),(1,1,-1),(1,-1,1),(1,1,1)]
    colors=[(0,0,0),(200,0,0),(0,200,0),(0,0,200), (200,200,0), (200,0,200), (0,200,200), (200,200,200)]
    #permutations=itertools.permutations(colors[1:])
    permutations=intelligentPermutations(colors)
    file_path="hex_chirality"
    distinct=[]
    doneConfigs=set()
    #pdb.set_trace()
    for j,perm in enumerate(permutations):
        #current_colors=[colors[0]]+list(perm)
        current_colors=perm
        rotations=rotate(current_colors)
        found=False
        for rotation in rotations:
            if rotation in doneConfigs:
                found=True
                break
        if not found:
            distinct.append(current_colors)
            for rotation in rotations:
                doneConfigs.add(rotation)
    for j,current_colors in enumerate(distinct):
        pointList=[]
        faceList=[]
        for i,p in enumerate(points):
            point(p[0],p[1],p[2], current_colors[i], pointList)
        face_points=[(0,2,3,1),(2,6,7,3),(6,4,5,7),(4,0,1,5),(0,4,6,2),(3,7,5,1)]
        for fp in face_points:
            face4(fp[0], fp[1], fp[2], fp[3], faceList)
        outfile=open(file_path+str(j)+".ply", "w")
        outfile.write("ply\n")
        outfile.write("format ascii 1.0\n")
        outfile.write("obj_info hex\n")
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

def face4(a,b,c,d,facelist):
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

def face3(a,b,c,facelist):
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
        
def point(x,y,z,color,pointlist):
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

#tetChirality()
hexChirality()