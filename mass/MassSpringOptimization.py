import numpy as np
import networkx as nx
import pylab as plt
import time
import pdb
import math

radius=6.0
radius2=2*radius
radius4=4*radius
radius3=radius**3
atan32=math.atan(-3.2)

def handle_spring_constraints(network,cellList, jointsNetwork, pendingJoints, jointLookup, noDetach, detachedJoints):
    """
    """
    forces=[]
    for i in range(len(cellList)):
        forces.append([0.0,0.0,0.0])
    #if len(cellList)==2:
    #    pdb.set_trace()
    k_pull=.0000015
    k_push=.0000032
    global radius
    global radius2
    global radius4
    global radius3
    global atan32
    col = 0 #active collision/overlap error
    optimal_distance = 11.0 #slighly less than 2*radius. dont want to overcorrect and loose joints
    N_overlapping_pairs = 0
    compression=dict()
    #edges = network.edges()
    #for edge in network.edges():
    for edge in network:
        index_1 =  edge[0]
        index_2 = edge[1]
        obj1 = cellList[index_1]
        obj2 = cellList[index_2]

        if(index_1 != index_2):
            #that is to say these are NOT internal constraints pciked up
            #by the nearest neighbor approach
            v12 = obj2["location"] - obj1["location"]
            dist = Magnitude(v12)
            unit_norm_vector = NormVec(v12)
            cell_overlap = (optimal_distance-dist)
            id1=obj1["id"]
            id2=obj2["id"]
            if jointsNetwork[id1][id2]:
                if id1<id2:
                    pendingJoints[(id1,id2)]=jointLookup[id1][id2][1]
                else:
                    pendingJoints[(id2,id1)]=jointLookup[id1][id2][1]
                joint=jointLookup[id1][id2][0]
                joint.detach()
            #check for a collision
            if(dist < optimal_distance):
                #then apply the collision
                if id1 not in compression:
                    compression[id1]=[1.0,[]]
                if id2 not in compression:
                    compression[id2]=[1.0,[]]
                shared_volume=.0625*((radius4+dist)*(radius2-dist)**2/radius3-(radius4+optimal_distance)*(radius2-optimal_distance)**2/radius3)
                compression[id1][0]+=-shared_volume
                compression[id2][0]+=-shared_volume
                compression[id1][1].append([id2,unit_norm_vector])
                compression[id2][1].append([id1,-unit_norm_vector])
                #report updated distance
                #v12 = (obj2["displacement"]+obj2["location"]) - (obj1["displacement"]+obj1["location"])
                #dist_after = Magnitude(v12)
                #print("(%i-%i) Dist before:%.02f, Dist after: %.02f:" %(index_1,index_2,dist,dist_after))
                col += cell_overlap
                N_overlapping_pairs+=1
            elif (dist > optimal_distance):
                # apply attractive force. Rervse of ruplision vector d
                #d = -unit_norm_vector*(cell_overlap/2.0)
                #td changed to:
                d = unit_norm_vector*(cell_overlap/2.0)
                if dist>radius*(2.2+0.85*(math.atan(6.0*jointLookup[id1][id2][1]-3.2)-atan32)):
                    #if id2 not in jointLookup[id1]:
                    #    pdb.set_trace()
                    if id1<id2:
                        if (id1,id2) not in noDetach:
                            jointsNetwork[id1][id2]=False
                            jointsNetwork[id2][id1]=False
                            joint=jointLookup[id1].pop(id2)
                            jointLookup[id2].pop(id1)
                            if (id1,id2) in pendingJoints:
                                pendingJoints.pop((id1,id2))
                            else:
                                joint.detach()
                    else:
                        if (id2,id1) not in noDetach:
                            jointsNetwork[id1][id2]=False
                            jointsNetwork[id2][id1]=False
                            joint=jointLookup[id1].pop(id2)
                            jointLookup[id2].pop(id1)
                            if (id2,id1) in pendingJoints:
                                pendingJoints.pop((id2,id1))
                            else:
                                joint.detach()
                    detachedJoints.append((id1,id2))
                else:
                    #add this to the collision vec
                    #obj1["displacement"] = obj1["displacement"] - d*k_pull
                    #obj2["displacement"] = obj2["displacement"] + d*k_pull
                    
                    #switch to force model instead of displacement
                    force_vector=d*k_pull*jointLookup[id1][id2][1]
                    #print(force_vector, v12)
                    for i in range(3):
                        forces[id1][i]=forces[id1][i]-force_vector[i]
                        forces[id2][i]=forces[id2][i]+force_vector[i]
                    #if (force_vector[0]**2+force_vector[1]**2+force_vector[2]**2)**.5>9e-7:
                    #    pdb.set_trace()
                    #report updated distance
                    #v12 = (obj2["displacement"]+obj2["location"]) - (obj1["displacement"]+obj1["location"])
                    #dist_after = Magnitude(v12)
                    #print("(%i-%i) Dist before:%.02f, Dist after: %.02f:" %(index_1,index_2,dist,dist_after))
                    #col += abs(cell_overlap)
                    N_overlapping_pairs+=1
    for id1 in compression:
        vol_frac=compression[id1][0]
        #set a max force to deal with issues where fraction is found to be negative (because of overlapping impinging cells)
        #and to prevent forces from getting too large.
        if vol_frac<.001:
            vol_frac=.001
        neighbors=compression[id1][1]
        for neighbor in neighbors:
            force_vector=(k_push/vol_frac-k_push)*neighbor[1]
            #if (force_vector[0]**2+force_vector[1]**2+force_vector[2]**2)**.5>9e-7:
            #    pdb.set_trace()
            #vectorToNeighbor=[cellList[neighbor[0]]["location"][0]-cellList[cell]["location"][0], 
                              #cellList[neighbor[0]]["location"][1]-cellList[cell]["location"][1],
                              #cellList[neighbor[0]]["location"][2]-cellList[cell]["location"][2]]
            #print(vectorToNeighbor, force_vector)
            id2=neighbor[0]
            for i in range(3):
                forces[id2][i]=forces[id2][i]+force_vector[i]
    return forces
            #cellList[neighbor[0]]["displacement"]=cellList[neighbor[0]]["displacement"]+(k_push/vol_frac-k_push)*neighbor[1]
    #Apply displacemnts
    #for cell in cellList:
        #print("Displacement for cell %i is"% cell["id"],cell["displacement"])
        #cell["location"] = cell["location"] + cell["displacement"]
        #cell["displacement"] = np.array([0,0,0])

    #network.remove_edges_from(remove_edges)
    #print("Total Overlap %.03f\n%i overlapping parirs" %(col,N_overlapping_pairs))
    #if N_overlapping_pairs>0:
    #    avg_collision_error = col / (N_overlapping_pairs)
    #else:
        #avg_collision_error = 0
    #return avg_collision_error

def ScaleVec(v1, s):
    """ Scales a vector f*[x,y,z] = [fx, fy, fz]
        Returns - a new scaled vector [x,y,z] as a numpy array
    """
    return np.array([x*s for x in v1])

def Magnitude(v1):
    """ Computes the Magnitudenitude of a vector
        Returns - a float representing the vector magnitude
    """
    return np.sqrt(np.sum([x**2 for x in v1]))

def NormVec(v1):
    """ Computes a normalized version of the vector v1 or a unit vector.
        Returns - a normalizerd vector [x,y,z] as a numpy array
    """
    mag = Magnitude(v1)
    if(mag == 0):
        return np.array([0,0,0])
    return np.array([ x/mag  for x in v1 ])

def compute_interactionsOld(cellList, jointsNetwork):
    """ 
    Add edges between bordering cells in the graph.
    This updates cell neighbors between optimization steps.
    Note: interactions/joins occur at dist=12.0um, but collision/repulsion pushes
        cells away unitl approx 11.0um. This gives us room to add joints.
    """
    G = nx.Graph()
    interaction_length = 11.0

    for i in range(len(cellList)):
        for j in range(i+1,len(cellList)):
            distance = Magnitude(cellList[i]["location"]- cellList[j]["location"])
            if distance < interaction_length or jointsNetwork[cellList[i]["id"]][cellList[j]["id"]]:
                G.add_edge(i,j)
            

    #print("New edges are ",G.edges())
    return G

def compute_interactionsOld2(cellList, jointsNetwork):
    """ 
    Add edges between bordering cells in the graph.
    This updates cell neighbors between optimization steps.
    Note: interactions/joins occur at dist=12.0um, but collision/repulsion pushes
        cells away unitl approx 11.0um. This gives us room to add joints.
    """
    G = []
    interaction_length = 11.0

    for i in range(len(cellList)):
        for j in range(i+1,len(cellList)):
            distance = Magnitude(cellList[i]["location"]- cellList[j]["location"])
            if distance < interaction_length or jointsNetwork[cellList[i]["id"]][cellList[j]["id"]]:
                G.append((i,j))
            

    #print("New edges are ",G.edges())
    return G

def compute_interactions(cellList, jointLookup):
    edges=set()
    grid=dict()
    interaction_length = 11
    for i in range(len(cellList)):
        location=cellList[i]["location"]
        x=int(location[0]//interaction_length)
        y=int(location[1]//interaction_length)
        z=int(location[2]//interaction_length)
        if x not in grid:
            grid[x]=dict()
            grid[x][y]=dict()
            grid[x][y][z]=[i]
        elif y not in grid[x]:
            grid[x][y]=dict()
            grid[x][y][z]=[i]
        elif z not in grid[x][y]:
            grid[x][y][z]=[i]
        else:
            grid[x][y][z].append(i)
    for x in grid:
        for y in grid[x]:
            for z in grid[x][y]:
                if x+1 in grid:
                    for y2 in range(y-1, y+2):
                        if y2 in grid[x+1]:
                            for z2 in range(z-1, z+2):
                                if z2 in grid[x+1][y2]:
                                    for i in grid[x][y][z]:
                                        for j in grid[x+1][y2][z2]:
                                            distance = Magnitude(cellList[i]["location"]- cellList[j]["location"])
                                            if distance<interaction_length:
                                                edges.add((i,j))
                if y+1 in grid[x]:
                    for z2 in range(z-1, z+2):
                        if z2 in grid[x][y+1]:
                            for i in grid[x][y][z]:
                                for j in grid[x][y+1][z2]:
                                    distance = Magnitude(cellList[i]["location"]- cellList[j]["location"])
                                    if distance<interaction_length:
                                        edges.add((i,j))
                if z+1 in grid[x][y]:
                    for i in grid[x][y][z]:
                        for j in grid[x][y][z+1]:
                            distance = Magnitude(cellList[i]["location"]- cellList[j]["location"])
                            if distance<interaction_length:
                                edges.add((i,j))
                for a in range(len(grid[x][y][z])):
                    i=grid[x][y][z][a]
                    for b in range(a+1, len(grid[x][y][z])):
                        j=grid[x][y][z][b]
                        distance = Magnitude(cellList[i]["location"]- cellList[j]["location"])
                        if distance<interaction_length:
                            edges.add((i,j))
    for i in jointLookup:
        for j in jointLookup[i]:
            if (i,j) not in edges and (j,i) not in edges:
                edges.add((i,j))
    return edges
        


def OptimizeCellPositions(cellList, jointsNetwork, times, pendingJoints, jointLookup, noDetach, detachedJoints):
    """
    Function: use mass-spring system to reposition cells are mitosis.
    # Cannot use elastic collisions since these cells have no apparent velocity.
    """
    #print("In")
    #for x in cellList:
    #    print(x[0].getPosition())
    #print("Optimizing cell positions with Mass-Spring Optimization")
    cellList_local = []
    lasttime=time.time()
    # Convert ODE object into a python object
    for i in range(len(cellList)):
        cell_body,cell_geom = cellList[i]
        data = cell_body.getData()
        cell_id = data["id"]
        assert cell_id == i, "Cell ID does not match cell index. This is an assumption in this mass spting optimization"

        position = cell_body.getPosition()
        optimizaion_object = {"id": cell_id, "body":cell_body,
                              "location":np.array([position[0],position[1],position[2]]),
                              "displacement":np.array([0,0,0])}

        cellList_local.append(optimizaion_object)
    times[0]+=time.time()-lasttime
    lasttime=time.time()
    """
    a1 = {"id":0,"location":np.array([0,0]),
      "displacement":np.array([0,0])}

    a2 = {"id":1,"location":np.array([0,11]),
      "displacement":np.array([0,0])}

    a3 = {"id":2,"location":np.array([0,5.5]),
      "displacement":np.array([0,0])}

    cellList = [a1,a2,a3]
    #G = nx.Graph()
    #G.add_edges_from([(0,1),(0,2),(1,2)])
    """
    #do opitimizaiton
    avg_error_limit = 0.1
    #Tristan: set to 1 iteration because we will run more frequently
    max_iterations = 1
    iteration = 1
    avg_collision_error = 100.0

    #while avg_collision_error >= avg_error_limit and iteration <= max_iterations:
        # compute interaction network. quadtree or double for loop
    network = compute_interactions(cellList_local, jointLookup)
    times[1]+=time.time()-lasttime
    lasttime=time.time()
    # compute displacements from cell collisions.
    forces=handle_spring_constraints(network,cellList_local, jointsNetwork, pendingJoints, jointLookup, noDetach, detachedJoints)
    times[2]+=time.time()-lasttime
    lasttime=time.time()
    return forces
    #print("Average Overlap Error(i=%i):" % iteration,avg_collision_error)
    #print(cellList_local)
    

    """
    plt.figure()
    x_vec = [x["location"][0] for x in cellList]
    y_vec = [x["location"][1] for x in cellList]
    plt.scatter(x_vec,y_vec)
    plt.xlim([-20,20])
    plt.ylim([-20,20])
    plt.savefig("iteration%i.png" %iteration)
    """
    #iteration+=1

    # adjust positions in origional cell list oobject
    #for object_master,object_local in zip(cellList,cellList_local):
    #    new_pos = object_local["location"] #x,y,z as numpy array
    
        #cell_body, cell_geom = object_master
        #cell_body.setPosition(new_pos[0],new_pos[1],new_pos[2])
    #times[3]+=time.time()-lasttime
    #lasttime=time.time()
    #print("Out")
    #for x in cellList:
    #    print(x[0].getPosition())
    #return cellList