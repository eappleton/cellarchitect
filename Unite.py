import sys, pdb
sys.path.append("mesh")
sys.path.append("block_definition")
sys.path.append("developmental_tree")
sys.path.append("circuit_design")
sys.path.append("mass")
import Meshing
import Algorithms_for_MLST as mlst
import More_connection_algs as conAlgs
import NewDevTree as ndt
import DesignCircuits as circuits
import StateMachine
#import StateMachineDrawings as drawings
import parallelizeMass as sim
from decimal import Decimal as D


def mesh(meshType, shapeFile, maxCells, name, numCores):
    if meshType=="Hex":
        meshing=Meshing.HexMeshing()
    else:
        meshing=Meshing.TetMeshing()
    return meshing.run(shapeFile, name, maxCells, numCores)

def connections(adjacency, algorithm, numFaces):
    graph_mst= mlst.mlst(adjacency, numFaces)
    if graph_mst==dict():
        graph_mst=conAlgs.minimumSpanningTree(adjacency, numFaces)
    if algorithm == 'MLST':
        return graph_mst
    elif algorithm == 'freely_moving_part':
        return conAlgs.freeMovingParts2(graph_mst, adjacency)
    else:
        raise ValueError("Algorithm Type")

def devTree(meshType, adjacency, colorBlocks=False):
    if meshType=="Hex":
        return ndt.makeDevTree("Hex", adjacency, 6, colorBlocks)
    else:
        return ndt.makeDevTree("Tet", adjacency, 4, colorBlocks)

def makeCircuits(tree, numDivisions, counterArchitecture):
    register, counterRecs=circuits.designRegister(tree, numDivisions, counterArchitecture)
    counter=circuits.designCounter(numDivisions, counterRecs, counterArchitecture)
    return counter, register

def run(name, numCores, maxCells, shapeFile, meshType, colorBlocks, numSims=0):
    try:
        fullgraph=dict()
        #fullgraph[0]=[-1,1,11,-1,-1,-1]
        #fullgraph[1]=[0,2,-1,-1,-1,-1]
        #fullgraph[2]=[1,3,-1,-1,-1,-1]
        #fullgraph[3]=[2,-1,4,-1,-1,-1]
        #fullgraph[4]=[-1,-1,5,3,-1,-1]
        #fullgraph[5]=[-1,-1,6,4,-1,-1]
        #fullgraph[6]=[7,-1,-1,5,-1,-1]
        #fullgraph[7]=[8,6,-1,-1,-1,-1]
        #fullgraph[8]=[9,7,-1,-1,-1,-1]
        #fullgraph[9]=[-1,8,-1,10,-1,-1]
        #fullgraph[10]=[-1,-1,9,11,-1,-1]
        #fullgraph[11]=[-1,-1,10,0,-1,-1]
        
        #fullgraph[0]=[-1,1,7,-1,-1,-1]
        #fullgraph[1]=[0,2,-1,-1,-1,-1]
        #fullgraph[2]=[1,-1,3,-1,-1,-1]
        #fullgraph[3]=[-1,-1,4,2,-1,-1]
        #fullgraph[4]=[5,-1,-1,3,-1,-1]
        #fullgraph[5]=[6,4,-1,-1,-1,-1]
        #fullgraph[6]=[-1,5,-1,7,-1,-1]
        #fullgraph[7]=[-1,-1,6,0,-1,-1]
        
        #fullgraph[0]=[-1,1,-1,-1,-1,-1]
        #fullgraph[1]=[0,2,-1,-1,-1,-1]
        #fullgraph[2]=[1,-1,-1,-1,-1,-1]
        
        #fullgraph[1]=[-1,-1,-1,2,-1,-1]
        #fullgraph[2]=[-1,3,1,-1,-1,-1]
        #fullgraph[3]=[2,-1,-1,-1,-1,-1]
        
        #fullgraph[0]=[1,-1,-1,-1]
        #fullgraph[1]=[0,-1,2,-1]
        #fullgraph[2]=[-1,-1,1,-1]
        #pdb.set_trace()
        fullgraph=mesh(meshType, shapeFile, maxCells, name, numCores)
        #graph=connections(fullgraph, "freely_moving_part", 6)
        #tree, numDivisions=devTree("Hex", graph)
        if meshType=="Hex":
            graph=connections(fullgraph, "freely_moving_part", 6)
        else:
            graph=connections(fullgraph, "freely_moving_part", 4)
        tree, numDivisions=devTree(meshType, graph, colorBlocks)
        counter, register=makeCircuits(tree, numDivisions, "linear")
        fileName="example_circuits/"+name+"_register.csv"
        file=open(fileName, "w")
        for x in register:
            file.write(x+",\n")
        file.close()
        #drawings.diagramFromFile(fileName, fileName[:-3]+".png")
        fileName2="example_circuits/"+name+"_counter.csv"
        file=open(fileName2, "w")
        for x in counter:
            file.write(x+",\n")
        file.close()
        #drawings.diagramFromFile(fileName2, fileName2[:-3]+".png")
        if numSims>0:
            machine=StateMachine.StateMachine(["Simple"], "example_circuits/"+name+"_state_machine_002.txt", fileName2, fileName, "High")
            machine.run(.002, numCores)
            p=sim.Parallelize()
            p.run(name+"_mass_simulation", "example_circuits/"+name+"_state_machine_002_simple.csv", numSims, numCores)
    except Meshing.NotConnectedError:
        print("Unable to mesh the shape with so few cells.")
        
def runExamples(numCores, meshType, numSims=0, numDays=7, sim_save_frequency=D('60.0')):
    examples=[]
    fullgraph=dict()
    fullgraph[0]=[-1,-1,-1,-1]
    #examples.append(["1Tet",fullgraph])
    fullgraph=dict()
    fullgraph[0]=[1,-1,-1,-1]
    fullgraph[1]=[0,-1,-1,-1]
    #examples.append(["2Tets", fullgraph])
    fullgraph=dict()
    fullgraph[0]=[1,-1,-1,-1]
    fullgraph[1]=[2,0,-1,-1]
    fullgraph[2]=[-1,-1,1,-1]
    #examples.append(["3Tets", fullgraph])
    fullgraph=dict()
    fullgraph[0]=[1,2,3,4]
    fullgraph[1]=[0,-1,-1,-1]
    fullgraph[2]=[0,-1,-1,-1]
    fullgraph[3]=[0,-1,-1,-1]
    fullgraph[4]=[0,-1,-1,-1]
    examples.append(["5Tets", fullgraph])
    fullgraph=dict()
    fullgraph[0]=[1,-1,-1,-1]
    fullgraph[1]=[0,-1,2,-1]
    fullgraph[2]=[-1,3,1,-1]
    fullgraph[3]=[-1,4,-1,2]
    fullgraph[4]=[-1,-1,-1,3]
    examples.append(["5Tets_line", fullgraph])
    for example in examples:
        name=example[0]
        fullgraph=example[1]
        #pdb.set_trace()
        #graph=connections(fullgraph, "freely_moving_part", 6)
        #tree, numDivisions=devTree("Hex", graph)
        if meshType=="Hex":
            graph=connections(fullgraph, "freely_moving_part", 6)
        else:
            graph=connections(fullgraph, "freely_moving_part", 4)
        tree, numDivisions, block_reference=devTree(meshType, graph, True)
        counter, register=makeCircuits(tree, numDivisions, "linear")
        fileName="example_circuits/"+name+"_register.csv"
        file=open(fileName, "w")
        for x in register:
            file.write(x+",\n")
        file.close()
        #drawings.diagramFromFile(fileName, fileName[:-3]+".png")
        fileName2="example_circuits/"+name+"_counter.csv"
        file=open(fileName2, "w")
        for x in counter:
            file.write(x+",\n")
        file.close()
        #drawings.diagramFromFile(fileName2, fileName2[:-3]+".png")
        if numSims>0:
            #machine=StateMachine.StateMachine(["Simple"], "example_circuits/"+name+"_state_machine_high_002.txt", fileName2, fileName, "High")
            #machine.run(.002, numCores)
            #p=sim.Parallelize()
            #p.run(name+"high_mass_simulation", "example_circuits/"+name+"_state_machine_high_002_simple.csv", numSims, numCores, numDays, block_reference, sim_save_frequency)
            machine=StateMachine.StateMachine(["Simple"], "example_circuits/"+name+"_state_machine_max_002.txt", fileName2, fileName, "Max")
            machine.run(.002, numCores)
            p=sim.Parallelize()
            p.run(name+"max_mass_simulation", "example_circuits/"+name+"_state_machine_max_002_simple.csv", numSims, numCores, numDays, block_reference, sim_save_frequency)
    

    
if __name__=="__main__":
    if 0==1:
        run("test_binary", 4, 20)
    if 0==1:
        run("3tet_high", 8, 16)
    if 0==1:
        run("SimpleFrameFullTest", 4, 512,"mesh/examples/simpleframe.STL","Hex", False)
    if 0==1:
        runExamples(12, "Tet", 96)
    if 1==1:
        