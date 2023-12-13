from multiprocessing import Process, Queue
import multiprocessing
import queue
import os
import random
from decimal import Decimal as D

#local module
from run_simulation import Run_Simulation
class Parallelize():
    def __init__(self):
        self.identifier=str(random.randint(1,1000000))
    
    def run(self, jobTitle, stateFile, numSims, numCores, numDays, block_reference, sim_save_frequency=D('15.0')):
        self.numDays=numDays
        self.save_frequency=sim_save_frequency
        self.block_reference=block_reference
        nextJob=Queue()
        processes=[]
        print("Starting "+str(numSims)+" simulations on "+str(numCores)+" cores.")
        for i in range(numSims):
            nextJob.put((jobTitle, stateFile, Run_Simulation()))
        print(nextJob.qsize())
        for i in range(numCores):
            processes.append(Process(target=self.massSim,args=(nextJob, str(i))))
            processes[i].start()
            print(id(processes[i]))
            print(processes[i])
        for p in processes:
            p.join()
            
    def runOnFolder(self, folderPath, numSims, numCores, numDays):
        self.numDays=numDays
        fileList=os.listdir(folderPath)
        nextJob=Queue()
        processes=[]
        for file in fileList:
            for i in range(numSims):
                nextJob.put((str(file)[:-4], folderPath+"/"+str(file), Run_Simulation()))
        for i in range(numCores):
            processes.append(Process(target=self.massSim, args=(nextJob, str(i))))
            processes[i].start()
            
    def massSim(self, nextJob, processID):
        while True:
            try:
                job=nextJob.get(block=True, timeout=1)
                sim=job[2]
                print("sim ",sim, id(sim))
                sim.begin(self.block_reference, job[0]+"_"+self.identifier+"_"+processID,job[1], self.numDays, self.save_frequency)
            except queue.Empty:
                print("Exit")
                return
    def makeSimulator(self):
        return Run_Simulation()
if __name__=='__main__':
    multiprocessing.set_start_method('fork')
    p=Parallelize()
    #p.runOnFolder("examples/Random/Batch1",10,20)
    #p.run("cone_hetero", "examples/stateMachine_cone_002_max_simple.csv", 8, 4)
    p.run("2019_02_20_3Cubes", "examples/test_rod_max_state_machine_002_simple.csv", 48, 10)