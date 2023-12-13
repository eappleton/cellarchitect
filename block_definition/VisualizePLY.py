
# coding: utf-8

# In[45]:


import random


# In[46]:


def cubePoints(x,y,z,color, pointlist):
    #top
    point(x-1.5,y+1.5,z+1.5, color, pointlist)
    point(x+1.5,y+1.5,z+1.5, color, pointlist)
    point(x+1.5,y+1.5,z-1.5, color, pointlist)
    point(x-1.5,y+1.5,z-1.5, color, pointlist)
    
    #bottom
    point(x-1.5,y-1.5,z+1.5, color, pointlist)
    point(x+1.5,y-1.5,z+1.5, color, pointlist)
    point(x+1.5,y-1.5,z-1.5, color, pointlist)
    point(x-1.5,y-1.5,z-1.5, color, pointlist)


# In[47]:


def point(x,y,z,color,pointlist):
    pointlist.append(str(x)+" "+str(y)+" "+str(z)+" "+str(color[0])+" "+str(color[1])+" "+str(color[2])+" \n")


# In[48]:


def face4(a,b,c,d,facelist):
    facelist.append("4 "+str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" \n")


# In[49]:


def cover(start, direction, facelist):
    if direction==0:
        face4(start+7,start+4, start, start+3,facelist)
    elif direction==1:
        face4(start+5,start+6, start+2, start+1,facelist)
    elif direction==2:
        face4(start+7,start+6, start+5, start+4,facelist)
    elif direction==3:
        face4(start,start+1, start+2, start+3,facelist)
    elif direction==4:
        face4(start+7,start+6, start+2, start+3,facelist)
    else:
        face4(start+4,start+5, start+1, start,facelist)


# In[50]:


def connectorParent(x, y, z, color, start, count, direction, pointlist, facelist):
    if direction==0:
        point(x-1.5,y-.5,z-.5, color, pointlist)
        point(x-1.5,y-.5,z+.5, color, pointlist)
        point(x-1.5,y+.5,z+.5, color, pointlist)
        point(x-1.5,y+.5,z-.5, color, pointlist)
        face4(start+4, start+7, count, count+1, facelist)
        face4(start+7, start+3, count+3, count, facelist)
        face4(start+3, start, count+2, count+3, facelist)
        face4(start, start+4, count+1, count+2, facelist)
    elif direction==1:
        point(x+1.5,y-.5,z+.5, color, pointlist)
        point(x+1.5,y-.5,z-.5, color, pointlist)
        point(x+1.5,y+.5,z-.5, color, pointlist)
        point(x+1.5,y+.5,z+.5, color, pointlist)
        face4(start+6, start+5, count, count+1, facelist)
        face4(start+5, start+1, count+3, count, facelist)
        face4(start+1, start+2, count+2, count+3, facelist)
        face4(start+2, start+6, count+1, count+2, facelist)
    elif direction==2:
        point(x-.5,y-1.5,z-.5, color, pointlist)
        point(x+.5,y-1.5,z-.5, color, pointlist)
        point(x+.5,y-1.5,z+.5, color, pointlist)
        point(x-.5,y-1.5,z+.5, color, pointlist)
        face4(start+6, start+7, count, count+1, facelist)
        face4(start+7, start+4, count+3, count, facelist)
        face4(start+4, start+5, count+2, count+3, facelist)
        face4(start+5, start+6, count+1, count+2, facelist)
    elif direction==3:
        point(x-.5,y+1.5,z+.5, color, pointlist)
        point(x+.5,y+1.5,z+.5, color, pointlist)
        point(x+.5,y+1.5,z-.5, color, pointlist)
        point(x-.5,y+1.5,z-.5, color, pointlist)
        face4(start+1, start, count, count+1, facelist)
        face4(start, start+3, count+3, count, facelist)
        face4(start+3, start+2, count+2, count+3, facelist)
        face4(start+2, start+1, count+1, count+2, facelist)
    elif direction==4:
        point(x+.5,y-.5,z-1.5, color, pointlist)
        point(x-.5,y-.5,z-1.5, color, pointlist)
        point(x-.5,y+.5,z-1.5, color, pointlist)
        point(x+.5,y+.5,z-1.5, color, pointlist)
        face4(start+7, start+6, count, count+1, facelist)
        face4(start+6, start+2, count+3, count, facelist)
        face4(start+2, start+3, count+2, count+3, facelist)
        face4(start+3, start+7, count+1, count+2, facelist)
    else:
        point(x-.5,y-.5,z+1.5, color, pointlist)
        point(x+.5,y-.5,z+1.5, color, pointlist)
        point(x+.5,y+.5,z+1.5, color, pointlist)
        point(x-.5,y+.5,z+1.5, color, pointlist)
        face4(start+5, start+4, count, count+1, facelist)
        face4(start+4, start, count+3, count, facelist)
        face4(start, start+1, count+2, count+3, facelist)
        face4(start+1, start+5, count+1, count+2, facelist)


# In[51]:


def connectorChild(x, y, z, color, start, count, parentcount, direction, pointlist, facelist):
    connectorParent(x, y, z, color, start, count, direction, pointlist, facelist)
    if direction<2:
        face4(count,count+1, parentcount, parentcount+1, facelist)
        face4(count+1, count+2, parentcount+3, parentcount, facelist)
        face4(count+2, count+3, parentcount+2, parentcount+3, facelist)
        face4(count+3, count, parentcount+1, parentcount+2, facelist)
    elif direction<4:
        face4(count,count+1, parentcount+2, parentcount+3, facelist)
        face4(count+1, count+2, parentcount+1, parentcount+2, facelist)
        face4(count+2, count+3, parentcount, parentcount+1, facelist)
        face4(count+3, count, parentcount+3, parentcount, facelist)
    else:
        face4(count,count+1, parentcount, parentcount+1, facelist)
        face4(count+1, count+2, parentcount+3, parentcount, facelist)
        face4(count+2, count+3, parentcount+2, parentcount+3, facelist)
        face4(count+3, count, parentcount+1, parentcount+2, facelist)


# In[52]:


def makeShape(adjacencypath, groupspath, outfilepath, name, seed=None, numgroups=0):
    #read MLST
    adjacency=open(adjacencypath)
    outfile=open(outfilepath,"w")
    graph=dict()
    for line in adjacency:
        first=line.find(" ")
        second=first+4
        third=line.find(",",second+1)
        fourth=line.find(",",third+1)
        fifth=line.find(",",fourth+1)
        sixth=line.find(",",fifth+1)
        seventh=line.find(",",sixth+1)
        eighth=line.find("]",seventh+1)
        neighbors=[int(line[second:third]),int(line[third+1:fourth]),int(line[fourth+1:fifth]), int(line[fifth+1:sixth]),int(line[sixth+1:seventh]),int(line[seventh+1:eighth])]
        key=int(line[0:first])
        graph[key]=[]
        for n in neighbors:
            graph[key].append(n)
    adjacency.close()
    groupsfile=open(groupspath)
    colors=dict()
    count=0
    colorfrac=255/(numgroups-1)
    for line in groupsfile:
        second=line.find("[")
        if seed:
            color=[0,int(count*colorfrac),random.randint(0,255)]
            count+=1
        else:
            color=[random.randint(0,255),random.randint(0,255),random.randint(0,255)]
        while True:
            first=second
            second=line.find(",",first+1)
            if second==-1:
                second=line.find("]",first+1)
                colors[int(line[first+1:second])]=color
                break
            colors[int(line[first+1:second])]=color
    colors[seed]=[200,0,0]
    groupsfile.close()
    #first item in graph
    top=list(graph)[0]
    facelist=[]
    pointlist=[]
    count=build(graph, colors, [0,0,0], top, 0, pointlist, facelist)
    print(len(pointlist), count)
    outfile.write("ply\n")
    outfile.write("format ascii 1.0\n")
    outfile.write("obj_info "+name+"\n")
    outfile.write("element vertex "+str(count)+"\n")
    outfile.write("property float x\n")
    outfile.write("property float y\n")
    outfile.write("property float z\n")
    outfile.write("property uchar red\n")
    outfile.write("property uchar green\n")
    outfile.write("property uchar blue\n")
    outfile.write("element face "+str(len(facelist))+"\n")
    outfile.write("property list uchar int vertex_indices\n")
    outfile.write("end_header\n")
    for line in pointlist:
        outfile.write(line)
    for line in facelist:
        outfile.write(line)
    outfile.close()
    
def makeShape2(graph, outfilepath, name,  groupspath="", seed=None, numgroups=0):
    outfile=open(outfilepath,"w")
    if groupspath:
        colors=dict()
        count=0
        colorfrac=255/(numgroups-1)
        groupsfile=open(groupspath)
        for line in groupsfile:
            second=line.find("[")
            if seed:
                color=[0,int(count*colorfrac),random.randint(0,255)]
                count+=1
            else:
                color=[random.randint(0,255),random.randint(0,255),random.randint(0,255)]
            while True:
                first=second
                second=line.find(",",first+1)
                if second==-1:
                    second=line.find("]",first+1)
                    colors[int(line[first+1:second])]=color
                    break
                colors[int(line[first+1:second])]=color
        colors[seed]=[200,0,0]
        groupsfile.close()
    else:
        colors=dict()
        for node in graph:
            colors[node]=[random.randint(0,255),random.randint(0,255),random.randint(0,255)]
    #first item in graph
    top=list(graph)[0]
    facelist=[]
    pointlist=[]
    count=build(graph, colors, [0,0,0], top, 0, pointlist, facelist)
    print(len(pointlist), count)
    outfile.write("ply\n")
    outfile.write("format ascii 1.0\n")
    outfile.write("obj_info "+name+"\n")
    outfile.write("element vertex "+str(count)+"\n")
    outfile.write("property float x\n")
    outfile.write("property float y\n")
    outfile.write("property float z\n")
    outfile.write("property uchar red\n")
    outfile.write("property uchar green\n")
    outfile.write("property uchar blue\n")
    outfile.write("element face "+str(len(facelist))+"\n")
    outfile.write("property list uchar int vertex_indices\n")
    outfile.write("end_header\n")
    for line in pointlist:
        outfile.write(line)
    for line in facelist:
        outfile.write(line)
    outfile.close()


# In[53]:


def build(graph, colors, coordinates, vertex, start, pointlist, facelist, parent=-2, parentcount=0):
    cubePoints(coordinates[0], coordinates[1], coordinates[2], colors[vertex], pointlist)
    count=start+8
    todolist=[]
    for i in range(6):
        if graph[vertex][i]==-1:
            cover(start, i, facelist)
        elif graph[vertex][i]==parent:
            connectorChild(coordinates[0], coordinates[1], coordinates[2], colors[vertex], start, count, parentcount, i, pointlist, facelist)
            count+=4
        else:
            connectorParent(coordinates[0], coordinates[1], coordinates[2], colors[vertex], start, count, i, pointlist, facelist)
            child=[graph[vertex][i],count]
            if i==0:
                coordinates2=[coordinates[0]-8,coordinates[1],coordinates[2]]
            elif i==1:
                coordinates2=[coordinates[0]+8,coordinates[1],coordinates[2]]
            elif i==2:
                coordinates2=[coordinates[0],coordinates[1]-8,coordinates[2]]
            elif i==3:
                coordinates2=[coordinates[0],coordinates[1]+8,coordinates[2]]
            elif i==4:
                coordinates2=[coordinates[0],coordinates[1],coordinates[2]-8]
            else:
                coordinates2=[coordinates[0],coordinates[1],coordinates[2]+8]
            child.append(coordinates2)
            todolist.append(child)
            count+=4
    for child in todolist:
        count=build(graph, colors, child[2], child[0], count, pointlist, facelist, vertex, child[1])
    return count

if __name__=="__main__":
    # In[54]:
    
    
    if 0==1:
        groups="/Users/Tristan/Documents/ELM_Project/simpleframe_test_Groups.txt"
        mlst="/Users/Tristan/Documents/ELM_Project/simpleframe_test_MLST.txt"
        outfile="/Users/Tristan/Documents/ELM_Project/simpleframe_test_ply.ply"
        makeShape(mlst, groups, outfile, "simpleframe")
    
    
    # In[55]:
    
    
    #run on symmetry breakdown of tubesphere.
    if 0==1:
        groups="/Users/Tristan/Documents/ELM_Project/tubesphere_test_groups.txt"
        mlst="/Users/Tristan/Documents/ELM_Project/tubesphere_test_MLST.txt"
        outfile="/Users/Tristan/Documents/ELM_Project/tubesphere_test_ply.ply"
        makeShape(mlst, groups, outfile, "tubesphere", 1, 137)
    
    
    # In[56]:
    
    
    #run on symmetry breakdown of tubesphere.
    if 1==1:
        groups="/Users/Tristan/Documents/ELM_Project/tubespherebig_test_groups.txt"
        mlst="/Users/Tristan/Documents/ELM_Project/tubespherebig_test_MLST.txt"
        outfile="/Users/Tristan/Documents/ELM_Project/tubespherebig_test_ply.ply"
        makeShape(mlst, groups, outfile, "tubespherebig", 1, 35)
    
