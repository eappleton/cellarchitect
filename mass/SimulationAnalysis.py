"""
Helper functions for reporting at the end of simulations
"""
import os
import numpy as np
from scipy.spatial.distance import euclidean
import time

def report_locations(cellList):
    print("######################")
    print("Cell Position Reporter")
    print("######################")
    
    for idx,cell_object in enumerate(cellList):
        cell_body,cell_geom = cell_object
        position = cell_body.getPosition()
        
        data = cell_body.getData()
        cell_id = data["id"]
        cell_type = data["type"]

        # Report cell
        print("ID:%s, State:%s, Position:" % (cell_id,cell_type),position)
        # u,v,w = cell_body.getLinearVel()

def save_locations(cellList, jointsNetwork, expressedIgnored, block_reference, states, total_time,output_folder,fileprefix="cell_positions"):
    """
    Save cell positions in file as [file].csv.[time]
    Currently saves:
        x,y,z (position)
        cell_id (unique identifier for each cell)
        cell_state (identifier for a possible genetic makeup of a cell)
        radius
        color - color of rfps (or assigned based on cell state) expressed as 24bit rgb value converted to decimal integer.
        connections - cell_id the current cell is physically connected to.
    """
    lookup_time=0.0
    time1 = int(total_time)
    filename = "%s.csv.%s" % (fileprefix,str(time1))
    filepath = os.path.join(output_folder,filename)
    header = "x,y,z,cell_id,cell_state,radius,color,connections, block_id\n"
    colorlookup={"eGFP":(0,200,50), "mScarlet":(200,0,0), "iRFP670":(100,0,0), "tdTomato":(255,138,56), "EBFP2":(0,246,255), "EYFP":(221,255,73), "mBanana": (248,252,32), "CellDeath":(0,0,0), "mPlum":(226,66,244)}
    with open(filepath,'w') as outfile:
        # Write header
        outfile.write(header)
        bindingCounts=dict()
        totalBinding=0
        for idx,cell_object in enumerate(cellList):
            cell_body,cell_geom = cell_object
            position = cell_body.getPosition()
            
            data = cell_body.getData()
            cell_id = data["id"]
            cell_type = data["type"]
            radius = "6.0"
            joints=""
            neighbors=jointsNetwork[cell_id]
            for neighbor in neighbors:
                if neighbors[neighbor]:
                    joints=joints+str(neighbor)+";"
                    cell_type2=cellList[neighbor][0].getData()["type"]
                    totalBinding+=1
                    if cell_type<cell_type2:
                        if cell_type not in bindingCounts:
                            bindingCounts[cell_type]=dict()
                        if cell_type2 not in bindingCounts[cell_type]:
                            bindingCounts[cell_type][cell_type2]=1
                        else:
                            bindingCounts[cell_type][cell_type2]+=1
                    else:
                        if cell_type2 not in bindingCounts:
                            bindingCounts[cell_type2]=dict()
                        if cell_type not in bindingCounts[cell_type2]:
                            bindingCounts[cell_type2][cell_type]=1
                        else:
                            bindingCounts[cell_type2][cell_type]+=1
            red=0
            green=0
            blue=0
            count=0
            for gene in expressedIgnored[cell_type]:
                if gene in colorlookup:
                    red+=colorlookup[gene][0]
                    green+=colorlookup[gene][1]
                    blue+=colorlookup[gene][2]
                    count+=1
            if count>0:
                red=int(red/count)
                green=int(green/count)
                blue=int(blue/count)
                color=colorConvert(red, green, blue)
            else:
                color=colorConvert(150,150,150)
            proteins=states[cell_type]
            time3=time.time()
            current_dict=block_reference
            block_id="-1"
            #print(proteins)
            for protein in proteins:
                if current_dict and protein in current_dict:
                    block_id=str(current_dict[protein][0])
                    current_dict=current_dict[protein][1]
                #if the protein is also a block forming protein, but is not a match for the previous block forming proteins, then we set
                #block_id to -1 because the cell does not belong to a single block.
                elif protein in block_reference:
                    block_id="-1"
                    break
            lookup_time+=time.time()-time3
            line = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (str(position[0]),str(position[1]),str(position[2]),str(cell_id),str(cell_type),radius,color,joints, block_id)
            outfile.write(line)
        try:
            os.mkdir(output_folder+"/joint_stats/")
        except FileExistsError:
            pass
        outfile2=open(output_folder+"/joint_stats/joints%s.csv" % (str(time1)),"w")
        for key in bindingCounts:
            for key2 in bindingCounts[key]:
                outfile2.write(str(key)+","+str(key2)+","+str(bindingCounts[key][key2])+",\n")
        # Report cell
    #print(str(time.time())+"\tSaved cell positions to %s" % filename)
    return lookup_time

def colorConvert(red, green, blue):
    colorbits=""
    x=128
    while x>=1:
        if red>=x:
            red=red-x
            colorbits=colorbits+"1"
        else:
            colorbits=colorbits+"0"
        x=x/2
    x=128
    while x>=1:
        if green>=x:
            green=green-x
            colorbits=colorbits+"1"
        else:
            colorbits=colorbits+"0"
        x=x/2
    x=128
    while x>=1:
        if blue>=x:
            blue=blue-x
            colorbits=colorbits+"1"
        else:
            colorbits=colorbits+"0"
        x=x/2
    color=0
    powers2=[8388608,4194304,2097152,1048576,524288,262144,131072,65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1]
    for i,bit in enumerate(colorbits):
        color+=int(bit)*powers2[i]
    return str(color)

def report_distances(cellList):
    print("######################")
    print("Cell Distance Reporter")
    print("######################")
    
    for i in range(0,len(cellList)):
        for j in range(i+1,len(cellList)):
            # Cell 1 position
            cell_body1,cell_geom1 = cellList[i]
            position1 = cell_body1.getPosition()
            id1= cell_body1.getData()["id"]

            # Cell 2 position
            cell_body2,cell_geom2 = cellList[j]
            position2 = cell_body2.getPosition()
            id2= cell_body2.getData()["id"]

            # Euclidean Distance
            distance = euclidean(position1,position2)
            
            # Report cell
            print("ID1:%s,ID2:%s, Distance:" % (id1,id2),distance)
            # u,v,w = cell_body.getLinearVel()