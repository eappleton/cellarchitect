"""
Author: Demarcus Briers, Boston University
    BU Robotics Lab, Hybrid and NEtworked Systems Group
Function: 

Features:
    Performs manual collision detction.
    Calculates cell-cell adhesions.

TODO:
    boundary conditions are local to function
    boundry force is static and independent of dt
    jointsnetwork needs to have total cell count or more
    interaction network is manually created. needs to have chiose to read from FSM
"""
import pdb
import sys
#from direct.showbase.ShowBase import ShowBase
from panda3d.ode import OdeWorld,OdeBody, OdeMass, OdeSphereGeom
from panda3d.ode import OdeSimpleSpace, OdeHashSpace,OdeJointGroup,OdeBallJoint,OdeSliderJoint,OdeUtil,OdeContactJoint,OdeContact
from panda3d.core import BitMask32, CardMaker, Vec4, Quat

#import networkx as nx
from numpy.random import uniform
import numpy as np
from scipy.spatial.distance import euclidean
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal as D
from collections import defaultdict
import inspect, csv
import random, time
import math

#local modules
from MassSpringOptimization import *
import SimulationAnalysis
from FiniteStateMachine import FiniteStateMachine



class Run_Simulation():
    def begin(self, block_reference, PROJECT_LABEL=None, state_machine=None, days=11, save_frequency=D('15.0')):
        times={"collisions":0.0, "quickstep1": 0.0, "quickstep2": 0.0, "reform joints": 0.0, "standard movement": 0.0, "save": 0.0}
        deathmode=False
        self.optimizationTime=[0.0,0.0,0.0,0.0]
        self.noDetach=dict()
        self.readBindingStrengths()
        self.pi4=math.pi/4
        #pdb.set_trace()
        #base = ShowBase()
        np.random.seed()
        random.seed()
        ################################
        ### Simulation Settings
        ###############################
        if not PROJECT_LABEL:
            PROJECT_LABEL = sys.argv[1]
            state_machine = sys.argv[2]
        printCount=0
        TIME_STOP = 60*24*days #60 minutes * 20hours
        #TIME_STOP=60*12*13
        #TIME_STOP=60*24*2
        self.dt = D("0.10") #time step, minutes? (original comment said hrs but doesn't make sense)
        nextAdjustment=30
        # Physics world setup
        self.hasJoint = False
        distances = []
        self.num_collisions = 0
        self.jointLookup=dict()
        self.connected=dict()
        pendingJoints=dict()
        # Setup physics world for Open Dynamics Engine (ODE)
        self.odeWorld = OdeWorld()
        #print("odeworld ",id(self.odeWorld))
        self.odeWorld.setGravity(0, 0.0,0) #x,y,z
        self.odeWorld.initSurfaceTable(1) #num of surfaces
        self.odeWorld.setSurfaceEntry(0, 0, 150, 1.0, 0.2, 0.9, 0.00001, 0.0, 0.002)
        self.odeWorld.setErp(0.8)
        # Create a space and add a contactgroup to it to add the contact joints
        self.space = OdeHashSpace()
        contactgroup = OdeJointGroup()
        #WORLD_BOUNDARY = 916.0 # micrometers #temp changed to 1.5mm on 3/1/19
        WORLD_BOUNDARY =150 #4/1/19 changed back to 150 for sims for Noushin
        WORLD_BOUNDARY5=5*WORLD_BOUNDARY
        MAX_VELOCITY = 200.0 # micrometers/minute   
    
        #initialize cell populatiton
        self.cellList = []
        cell_interactions = nx.Graph()
        init_pos = [(25,0,0),(0,25,0)]
        init_f = [(0,0,0),(0,-2e-6,0)]
        for idx,cell_id in enumerate([0]):
            cell_obj = self.create_cell(self.odeWorld,init_pos[idx],init_f[idx],idx,0)
            self.cellList.append(cell_obj)
    
    
        ##########################################
        # Setup interactions and joints networks.
        #########################################3
        #FSM = FiniteStateMachine(file="stateMachine4count_2tetexc01_simple.csv")
        self.FSM = FiniteStateMachine(state_machine, self.bindingStrengths)
        OUTPUT_FOLDER = self.FSM.generate_sim_folder(state_machine=state_machine
                                                ,project_label=PROJECT_LABEL,
                                                output_folder="simulations_1.0")
        self.JointsNetwork = self.createJointNetwork(1536) #active cell adhesions due to collisions in space
        self.InteractionsNetwork = self.FSM.add_interaction_energies_fast(self.cellList) #allowed cell adhesions deermined by surface proteins
        
        #setup logic for first cell
        child_body = self.cellList[0][0]
        child_data = child_body.getData()
        new_child_data = self.FSM.updateAttributes(child_data,new_state_id=0)
        child_body.setData(new_child_data)
    
        ####################################
        # Run the simulation...
        ###################################
        current_time = D("0.0")
        nextSave=current_time+save_frequency
        #position_vec = [[],[]]
        count=0
        lookup_time=0.0
        while current_time<=TIME_STOP:
            #perform current time step
            current_time+=self.dt
    
            #########################################
            # perform mitosis. increment counter
            self.mitosis()
            #########################################
            time2=time.time()
            # Step the simulation and set the new positions
            self.odeWorld.quickStep(self.dt)
            if count%2==1:
                times["quickstep1"]=times["quickstep1"]+time.time()-time2
            else:
                times["quickstep2"]=times["quickstep2"]+time.time()-time2
            time1=time.time()
            self.space.collide((self.odeWorld,contactgroup), self.near_callback)
            time2=time.time()
            times["collisions"]=times["collisions"]+time2-time1
            
            
            #####################
            # BROWNIAN MOTION
            ###################
            if count%2==0:
                time1=time.time()
                for joint in pendingJoints:
                    j = OdeBallJoint(self.odeWorld)
                    j.attach(self.cellList[joint[0]][0], self.cellList[joint[1]][0])
                    self.jointLookup[joint[0]][joint[1]]=(j,pendingJoints[joint])
                    self.jointLookup[joint[1]][joint[0]]=(j,pendingJoints[joint])
                pendingJoints=dict()
                time2=time.time()
                times["reform joints"]=times["reform joints"]+time2-time1
                force_set=set()
                for i,cell_obj in enumerate(self.cellList):
                    cell_body,cell_geom = cell_obj
                    
                    #rendering
                    #cell_geom.setPosQuat(cell_body.getPosition(), Quat(cell_body.getQuaternion()))
                    x,y,z = cell_body.getPosition()
                    if abs(x)>WORLD_BOUNDARY5 or abs(y)>WORLD_BOUNDARY5 or abs(z)>WORLD_BOUNDARY5:
                    #    pdb.set_trace()
                         print("Extreme cell boundary violaion",x,y,z)
                         u,v,w = cell_body.getLinearVel()
                         print("Velocity at Extreme violation:",[u,v,w])
                         deathmode=True
                    #u,v,w = cell_body.getLinearVel()
                    #print(idx,":",u,v,w)
                    # reset veolocity
                    # Assume overdamped system. No accelertation
                    # position controlled by brownian motion.
                    #cell_body.setLinearVel(0.,0.,0.)
                    #cell_body.setForce(0.,0.,0.)
                    if i not in force_set:
                        # apply new random movemnt force
                        #relative to Demarcus's setting, halved force on 1/23/19 1:29pm
                        fx,fy,fz = self.random_walk_force(self.dt,1e-6)
                        #changed to somewhat consistent forces within shapes on 1/28/19 1pm
                        for j in self.connected[i]:
                            force_set.add(j)
                            cellbody2=self.cellList[i][0]
                            cellbody2.setForce(fx*(random.random()+.5),fy*(random.random()+.5),fz*(random.random()+.5))
                        #cell_body.setForce(fx,fy,fz) #used addForce(fx,fy,fz) in the past. possibly unsatble
                        cell_body.setForce(fx*(random.random()+.5),fy*(random.random()+.5),fz*(random.random()+.5))
            
                        # Apply boundary conditions.
                        # remove any brownian motion.
                    in_violation,vx,vy,vz = self.boundary_collision(cell_body,self.dt,WORLD_BOUNDARY)
                    if in_violation == True:
                        cell_body.setForce(0.,0.,0.)
                        cell_body.setLinearVel(vx,vy,vz)
                    
                    # Threshold velocities to prevent unstable and Zeno Behaviors
                    #max_velocity_violation,vx_2,vy_2,vz_2 = self.threshold_velocity(cell_body,MAX_VELOCITY)
                    #if max_velocity_violation == True:
                    #    cell_body.setForce(0.,0.,0.)
                    #    cell_body.setLinearVel(vx_2,vy_2,vz_2)
                self.applyDrag()
                times["standard movement"]=times["standard movement"]+time.time()-time2
            else:
                detachedJoints=[]
                forces=OptimizeCellPositions(self.cellList, self.JointsNetwork, self.optimizationTime, pendingJoints, self.jointLookup, self.noDetach, detachedJoints)
                for id1, id2 in detachedJoints:
                    self.fixConnected(id1,id2)
                for force in forces:
                    mag=Magnitude(force)
                    #adjustment for max position change of 2um
                    #if mag>3.619e-7:
                    #    for i in range(3):
                    #        force[i]=force[i]*(3.619e-7/mag)
                    #decreased to 1e-7 from 1.5e-7 on 2/20/19 3:12pm
                    if mag>1.0e-7:
                        for i in range(3):
                            force[i]=force[i]*(1.0e-7/mag)
                #print("Forces:",forces)
                for id1 in range(len(forces)):
                    self.cellList[id1][0].setForce(forces[id1][0],forces[id1][1],forces[id1][2])
            count+=1
            #if current_time>=nextAdjustment:
                #nextAdjustment+=30
                #OpitmizeCellPositions(cellList, self.JointsNetwork, self.optimizationTime)
            ###############
            #clean up
            ###############
            contactgroup.empty() # Clear the contact joints
    
            ##############################
            # Save positions at intervals for visualization.
            ################################
            time2=time.time()
            #time_stop_hours = TIME_STOP//60 + save_frequency # inlcude end item
            #if current_time==2340:
            #    pdb.set_trace()
            if deathmode:
                lookup_time+=SimulationAnalysis.save_locations(self.cellList, self.JointsNetwork, self.FSM.expressedIgnored, block_reference,self.FSM.states,
                                                  current_time*10,
                                                  output_folder=OUTPUT_FOLDER)
                print("lookup_time",lookup_time)
            else:
                if current_time>=nextSave:
                    nextSave+=save_frequency
                    if printCount==7:
                        print(count,self.optimizationTime)
                        print(times)
                        printCount=0
                    else:
                        printCount+=1
                    #print("Save time ")
                    #print(PROJECT_LABEL,"cell list", [(str(x[0].getData()), str(x[0].getPosition())) for x in cellList])
                    if current_time>=5000:
                        lookup_time+=SimulationAnalysis.save_locations(self.cellList, self.JointsNetwork, self.FSM.expressedIgnored,block_reference,self.FSM.states,
                                                      current_time*10,
                                                      output_folder=OUTPUT_FOLDER)
                        print("lookup_time",lookup_time)
                #elif 60*24*3.8>current_time>=60*24*3.3:
                elif current_time>=60*24*4.5:
                    lookup_time+=SimulationAnalysis.save_locations(self.cellList, self.JointsNetwork, self.FSM.expressedIgnored,block_reference,self.FSM.states,
                                                      current_time*10,
                                                      output_folder=OUTPUT_FOLDER)
                    #print("lookup_time",lookup_time)
            times["save"]+=time.time()-time2
            remove_noDetach=[]
            for idx in self.noDetach:
                self.noDetach[idx]=self.noDetach[idx]+1
                if self.noDetach[idx]==1:
                    remove_noDetach.append(idx)
            for idx in remove_noDetach:
                self.noDetach.pop(idx)
    
        #################################
        # End of simulation reporting
        #################################
        #print(min(distances))
        #t = range(len(distances))
        #plt.plot(t,distances)
        #plt.savefig('paired-distance.png',dpi=150)
    
        
        # plot single cell path
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        data = np.array(position_vec[0])
        ax.plot(data[:,0],data[:,1],data[:,2])
        ax.plot([data[0,0]],[data[0,1]],[data[0,2]], 'go')
        ax.plot([data[-1,0]], [data[-1,1]],[data[-1,2]], 'ro')
        plt.title("Cell 1 Path")
        plt.savefig("cell_1_path.png",dpi=150)
        """
        # report colliions
    
        # report locations
        #SimulationAnalysis.report_locations(cellList)
        #report distances
        #SimulationAnalysis.report_distances(cellList)
        print(self.optimizationTime)
        print(times)
    
    def readBindingStrengths(self):
        #checks to see if the constructor has been called from within this module or from Unite.py, because this changes the path
        #needed to find the recombinase efficiencies file.
        pyFile=inspect.stack()[-1].filename
        self.bindingStrengths=dict()
        try:
            if pyFile=="run_simulation.py" or pyFile=="parallelizeMass.py":
                reader=csv.reader(open("ProteinAffinities.csv"))
            else:
                reader=csv.reader(open("mass/ProteinAffinities.csv"))
            info=list(reader)
            #read in protein affinity Matrix
            protein_list=[]
            for protein in info[0][1:]:
                self.bindingStrengths[protein]=dict()
                protein_list.append(protein)
            for x in info[1:]:
                for i,y in enumerate(x[1:]):
                    self.bindingStrengths[x[0]][protein_list[i]]=float(y)
        except FileNotFoundError:
            print("No Protein Affinity File Found. Will Auto-Generate Affinities")
            return
    
    def create_cell(self, odeWorld,position,init_force,cell_id,cell_type):
        """
        Function: Create a biological cell.
        Requires a body(physics) and geometry(collisions)
        python tags/metadata ()
        """
        
        assert len(position) == 3, "Position vector for new cell is not 3D"
    
        # create spheres. 5 of them with radius 6.
        myBody = OdeBody(odeWorld)
        # Set mass. #REF:http://book.bionumbers.org/what-is-the-density-of-cells/
        myMass = OdeMass()
        myMass.setSphere(density=1e-12,radius=6.0) #units: grams/um^3 and um. true=1e-12
        myBody.setMass(myMass)
        myBody.setData({"id":cell_id,
                        "type":cell_type,
                        "mitosis_trigger":20*60, #20 hours, time is in minutes
                        "mitosis_counter":0})
    
        # Create a box geom for collision detection
        cell_geom = OdeSphereGeom(self.space, radius=6.0)
        cell_geom.setCollideBits(BitMask32.bit(0)) #collision masks
        cell_geom.setCategoryBits(BitMask32.bit(0)) #collision masks
        cell_geom.setBody(myBody)
    
        myBody.setPosition(position[0],position[1],position[2]) #x,y,z
        myBody.addForce(init_force[0],init_force[1],init_force[2]) #x,y,z
        self.connected[cell_id]={cell_id}
        return (myBody,cell_geom)
    
    def create_3Part_Joint(self, body1, body2, id1, id2):
        odeWorld=self.odeWorld
        pos1=body1.getPosition()
        pos2=body2.getPosition()
        vector=pos2-pos1
        mag=Magnitude(vector)
        if mag==0:
            return
        else:
            adjustment=mag/min(.1,mag/4.0)
            newPos1=vector/adjustment+pos1
            newPos2=vector/-adjustment+pos2
            newBody1=OdeBody(odeWorld)
            myMass = OdeMass()
            myMass.setSphere(density=1e-12,radius=1.0) #units: grams/um^3 and um. true=1e-12
            #cell_geom = OdeSphereGeom(self.space, radius=0.5)
            #cell_geom.setCollideBits(BitMask32.bit(0)) #collision masks
            #cell_geom.setCategoryBits(BitMask32.bit(0)) #collision masks
            #cell_geom.setBody(newBody1)
            newBody1.setMass(myMass)
            newBody1.setPosition(newPos1[0],newPos1[1],newPos1[2])
            #newBody1.setPosition(100,100,100)
            newBody2=OdeBody(odeWorld)
            myMass = OdeMass()
            myMass.setSphere(density=1e-12,radius=1.0) #units: grams/um^3 and um. true=1e-12
            #cell_geom = OdeSphereGeom(self.space, radius=0.5)
            #cell_geom.setCollideBits(BitMask32.bit(0)) #collision masks
            #cell_geom.setCategoryBits(BitMask32.bit(0)) #collision masks
            #cell_geom.setBody(newBody2)
            newBody2.setMass(myMass)
            newBody2.setPosition(newPos2[0],newPos2[1],newPos2[2])
            #newBody2.setPosition(200,200,200)
            j1=OdeBallJoint(odeWorld)
            j1.attach(body1,newBody1)
            j2=OdeBallJoint(odeWorld)
            j2.attach(body2, newBody2)
            j3=OdeSliderJoint(odeWorld)
            j3.attach(newBody1, newBody2)
            lowerID=min(id1,id2)
            higherID=max(id1,id2)
            if lowerID not in self.jointLookup:
                self.jointLookup[lowerID]=dict()
            self.jointLookup[lowerID][higherID]=(j1,j2,j3,newBody1,newBody2, body2)
        #pdb.set_trace()
    
    def applyDrag(self):
        drag_adjustment=5000000000.0
        for idx in range(len(self.cellList)):
            cell_body,cell_geom = self.cellList[idx] #cellList is a list of tuples.
            u,v,w = cell_body.getLinearVel()
            #for laminar flow, drag is proportional to velocity. Assumption of relatively laminar flow might be incorrect
            cell_body.addForce(-u/drag_adjustment, -v/drag_adjustment, -w/drag_adjustment)
    
    def near_callback(self, args, geom1, geom2):
        """Callback function for the collide() method.
    
        This function checks if the given geoms do collide and
        creates contact joints if they do.
        """
    
        # Check if the objects do collide
        #print("Collision occured!")
        #
        radius=6.0
        world,contactgroup = args
    
        body1=geom1.getBody()
        body2=geom2.getBody() 
        data1=body1.getData()
        data2=body2.getData()
        id1=data1["id"]
        id2=data2["id"]
        type1=data1["type"]
        type2=data2["type"]
        
        hasJoint = self.JointsNetwork[id1][id2]
        hasSurfaceBinding = False
        if self.InteractionsNetwork[type1][type2]>0:
            hasSurfaceBinding=True
    
        # create ballJoint if cells can adhere but havent came in contact yet.
        if hasJoint == False and hasSurfaceBinding == True:
            pos1=body1.getPosition()
            pos2=body2.getPosition()
            distance=((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)**.5
            canBind=False
            if distance<=2.2*radius:
                canBind=True
            #determine if cells are headed toward each other.
            else:
                displacement=pos2-pos1
                magd=Magnitude(displacement)
                v1=body1.getLinearVel()
                mag1=Magnitude(v1)
                if mag1!=0:
                    angle=math.acos(np.dot(displacement, v1)/(mag1*magd))
                if mag1==0 or -self.pi4<=angle<=self.pi4:
                    displacement=pos1-pos2
                    v2=body2.getLinearVel()
                    mag2=Magnitude(v2)
                    if mag2!=0:
                        angle=math.acos(np.dot(displacement, v2)/(mag2*magd))
                        if -self.pi4<=angle<=self.pi4:
                            canBind=True
                    elif mag1!=0:
                        canBind=True
                        
                            
            if canBind:
                self.num_collisions+=1
                j = OdeBallJoint(self.odeWorld)
                j.attach(body1, body2)
                if id1 not in self.jointLookup:
                    self.jointLookup[id1]=dict()
                if id2 not in self.jointLookup:
                    self.jointLookup[id2]=dict()
                self.jointLookup[id1][id2]=(j,self.InteractionsNetwork[type1][type2])
                self.jointLookup[id2][id1]=(j,self.InteractionsNetwork[type1][type2])
                if id1<id2:
                    self.noDetach[(id1,id2)]=0
                else:
                    self.noDetach[(id2,id1)]=0
                #inelastic collision. Only do if not already connected.
                if id2 not in self.connected[id1]:
                    count1=0
                    vx1,vy1,vz1=0.0,0.0,0.0
                    count2=0
                    vx2,vy2, vz2=0.0,0.0,0.0
                    for other in self.connected[id1]:
                        count1+=1
                        u,v,w=self.cellList[other][0].getLinearVel()
                        vx1+=u
                        vy1+=v
                        vz1+=w
                    for other in self.connected[id2]:
                        count2+=1
                        u,v,w=self.cellList[other][0].getLinearVel()
                        vx2+=u
                        vy2+=v
                        vz2+=w
                    vx_final=(vx1+vx2)/(count1+count2)
                    vy_final=(vy1+vy2)/(count1+count2)
                    vz_final=(vz1+vz2)/(count1+count2)
                    combined_connected=self.connected[id1].union(self.connected[id2])
                    for other in combined_connected:
                        self.connected[other]=combined_connected
                        self.cellList[other][0].setLinearVel(vx_final, vy_final, vz_final)
                #self.create_3Part_Joint(body1, body2, id1, id2)
                #j.setAnchor( (1,2,0) )
                #print("Adding joint",(id1,id2))
                self.JointsNetwork[id1][id2] = True
                self.JointsNetwork[id2][id1] = True
                #print("Adding joint")
        
        # compute elastic collision if cells dont bind.
        elif hasSurfaceBinding == False and hasJoint == False:
            self.num_collisions+=1
            #print("computing elastic collision") 
            self.elastic_collision(body1,body2, id1, id2)
        else:
            #print("Ignoring interaction") 
            pass
        
    
    def boundary_collision(self, cell_body,dt,boundary_size):
        """
        Function: Enforces fixed boundary conditions. Also thresholds cell velocity to prevent zeno behavior
        SEE: https://en.wikipedia.org/wiki/Hybrid_system
        """
        is_violation = False
        x,y,z = cell_body.getPosition()
        vx,vy,vz = cell_body.getLinearVel()
        x_adjust=self.check_boundary(x,vx,boundary_size)
        y_adjust=self.check_boundary(y,vy,boundary_size)
        z_adjust=self.check_boundary(z,vz,boundary_size)
        # Reverse cell velocities if outside bounding box
        if x_adjust==-1 or y_adjust==-1 or z_adjust==-1:
            is_violation = True
            vx = vx*x_adjust
            vy = vy*y_adjust
            vz = vz*z_adjust
    
    
        boundary_velocity = [is_violation,vx,vy,vz]
        return boundary_velocity
    
    def threshold_velocity(self, cell_body,MAX_VELOCITY):
        """
        # Threshold cell velocity. Prevent Zeno Behaviors when the velocity is very high.
        # if velocity is high cell can bounce between sides of the bounary infinitely.
        # SEE Zeno Behaviors in: https://en.wikipedia.org/wiki/Hybrid_system
        
        Returns 2 values:
            1) adjusted velocites
            2) a bool indicating if a cell exceeded the MAX_VELOCITY
        """
    
        vx,vy,vz = cell_body.getLinearVel()
        max_vel_violated = False
    
        if self.max_velocity_violated([vx,vy,vz],MAX_VELOCITY):
            vx = self.limit_velocity(vx,MAX_VELOCITY)
            vy = self.limit_velocity(vy,MAX_VELOCITY)
            vz = self.limit_velocity(vz,MAX_VELOCITY)
            max_vel_violated = True
    
    
        return [max_vel_violated,vx,vy,vz]
    
    def max_velocity_violated(self, velocities,MAX_VELOCITY):
        """ check if max velocity violated"""
        max_vel_violated = False
    
        for velocity in velocities:
            if velocity >= 0 and velocity > MAX_VELOCITY:
                max_vel_violated = True
            if velocity < 0 and velocity < -MAX_VELOCITY:
                max_vel_violated = True
    
        return max_vel_violated
    
    def limit_velocity(self, velocity,MAX_VELOCITY):
        """
        To prevent zeno behaviors, threshold max veolocity in the simulation.
    
        """
    
        if velocity >= 0:
            velocity = min(velocity,MAX_VELOCITY)
        elif velocity < 0:
             velocity = max(velocity,-MAX_VELOCITY)
    
        return velocity
    
    def check_boundary(self, x,vx,boundary_size):
        """
        Modify velocity if exceeeding boundaries.
        x - 1D poistion (float)
        vx - 1D velocity (float)
        boundary_size (int or float) min or max coordinate location in any direction. if 200
            then cell positions can range from -200,200 before boundary corrections are applied.
    
        Check if cell outside boundary and reverse velocity if so.
        Check is cell velocity has previously been reversed in the correct direction.
        Check if velocity above max velocity theshold.
        """
        #boundary_size = 200.0 #in one direction
        x_force = 0.0
    
        if x >= boundary_size and vx >= 0:
            x_force = -1.0
        elif x <= -boundary_size and vx <= 0:
            x_force = -1.0
        else:
            x_force = 1.0
    
        return x_force
    
    def random_walk_force(self, dt,scale=2e-6):
        """
        Gerate a 3D random walk movemnt.
        Uniform distribution is approximately brownian with small time steps.
        
        Inputs:
        ------------ 
        scale : magnitude of cell motlity for 1 unit of time (um/minute)
        dt - time step
    
        retur: 3 value array of x,y,z forces
        """
        
        #dt = .05
        #n = int(0/dt)
        #scale = 5.0
        random_walk = uniform(low=-scale,high=scale,size=(1,3))*float(dt)
        random_walk = random_walk[0]
    
        return random_walk
    
    def elastic_collision(self, body1,body2, id1, id2):
        """
        Compute velocity for body1 after an elastic collision with body2.
        
        REF:
        https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional_collision_with_two_moving_objects
        https://codereview.stackexchange.com/questions/125867/simulation-of-2d-elastic-balls
        
        TODO: consider the mass of aggregate objects. connected components.
        TODO: catch cases where norm is 0
        """
        if id2 not in self.connected[id1]:
            count1=0
            vx1,vy1,vz1=0.0,0.0,0.0
            count2=0
            vx2,vy2, vz2=0.0,0.0,0.0
            velos=dict()
            for other in self.connected[id1]:
                count1+=1
                velos[other]=self.cellList[other][0].getLinearVel()
                vx1+=velos[other][0]
                vy1+=velos[other][1]
                vz1+=velos[other][2]
            for other in self.connected[id2]:
                count2+=1
                velos[other]=self.cellList[other][0].getLinearVel()
                vx2+=velos[other][0]
                vy2+=velos[other][1]
                vz2+=velos[other][2]
            vx1=vx1/count1
            vy1=vy1/count1
            vz1=vz1/count1
            vx2=vx2/count2
            vy2=vy2/count2
            vz2=vz2/count2
            vx1f=((count1-count2)/(count1+count2))*vx1+(2*(count2)/(count1+count2))*vx2-vx1
            vy1f=((count1-count2)/(count1+count2))*vy1+(2*(count2)/(count1+count2))*vy2-vy1
            vz1f=((count1-count2)/(count1+count2))*vz1+(2*(count2)/(count1+count2))*vz2-vz1
            vx2f=((count2-count1)/(count1+count2))*vx2+(2*(count1)/(count1+count2))*vx1-vx2
            vy2f=((count2-count1)/(count1+count2))*vy2+(2*(count1)/(count1+count2))*vy1-vy2
            vz2f=((count2-count1)/(count1+count2))*vz2+(2*(count1)/(count1+count2))*vz1-vz2
            for other in self.connected[id1]:
                self.cellList[other][0].setLinearVel(vx1f+velos[other][0], vy1f+velos[other][1], vz1f+velos[other][2])
            for other in self.connected[id2]:
                self.cellList[other][0].setLinearVel(vx2f+velos[other][0], vy2f+velos[other][1], vz2f+velos[other][2])
        else:
            vx1, vy1, vz1= body1.getLinearVel()
            body1.setLinearVel(body2.getLinearVel())
            body2.setLinearVel(vx1,vy1,vz1)
    
    def cellInteractions(self, interaction=True):
        I = defaultdict(lambda:{})
        #for i in range(1):
        #    for j in range(4):
        #        I[str(i)][str(j)] = interaction
        #        I[str(j)][str(i)] = interaction
        I["0"]["0"] = True
        I["1"]["1"] = True
        
        I["1"]["0"] = False
        I["0"]["1"] = False
    
    
        return I
    
    def createJointNetwork(self, max_cells=200):
        G = defaultdict(lambda:{})
        for i in range(max_cells):
            for j in range(max_cells):
                if i>=j:
                    G[i][j] = False
                    G[j][i] = False
    
        return G
    
    def RandomPointOnSphere(self):
        """ Computes a random point on a sphere
            Returns - a point on a unit sphere [x,y,z] at the origin
        """
        u = np.random.random()*pow(-1., np.random.randint(0,1))
        theta = np.random.random()*2*np.pi
        x = np.sqrt(1-(u*u))*np.cos(theta)
        y = np.sqrt(1 - (u*u))*np.sin(theta)
        z = u
        return np.array((x,y,z))
    
    def fixConnected(self, id1, id2):
        set1=set()
        self.dfs(id1, set1)
        if id2 not in set1:
            self.connected[id1]=set1
            set2=set()
            self.dfs(id2, set2)
            self.connected[id2]=set2
            for other in self.connected[id1]:
                self.connected[other]=self.connected[id1]
            for other in self.connected[id2]:
                self.connected[other]=self.connected[id2]
        
    def dfs(self, id1, found):
        found.add(id1)
        for id2 in self.jointLookup[id1]:
            if id2 not in found:
                self.dfs(id2, found)
    
    def updateJoints(self, newCellStart, parentStates, parentIDs):
        radius=6
        oldInteractionsNetwork=self.InteractionsNetwork
        self.InteractionsNetwork=self.FSM.add_interaction_energies_fast(self.cellList)
        currentIndex=newCellStart
        for cell_body, cell_gemo in self.cellList[newCellStart:]:
            data1=cell_body.getData()
            id1=data1["id"]
            state1=data1["type"]
            if oldInteractionsNetwork[parentStates[parentIDs[id1]]][parentStates[parentIDs[id1]]]>0:
                cell_body2=self.cellList[parentIDs[id1]][0]
                data2=cell_body2.getData()
                state2=data2["type"]
                if self.InteractionsNetwork[state1][state2]>0:
                    id2=data2["id"]
                    j = OdeBallJoint(self.odeWorld)
                    j.attach(cell_body, cell_body2)
                    if id1 not in self.jointLookup:
                        self.jointLookup[id1]=dict()
                    if id2 not in self.jointLookup:
                        self.jointLookup[id2]=dict()
                    self.jointLookup[id1][id2]=(j, self.InteractionsNetwork[state1][state2])
                    self.jointLookup[id2][id1]=(j, self.InteractionsNetwork[state1][state2])
                    if id1<id2:
                        self.noDetach[(id1,id2)]=0
                    else:
                        self.noDetach[(id2,id1)]=0
                    #j.setAnchor( (1,2,0) )
                    #print("Adding joint",(id1,id2))
                    self.JointsNetwork[id1][id2] = True
                    self.JointsNetwork[id2][id1] = True
                    #print("Adding joint")
                    if id2 not in self.connected[id1]:
                        combined_connected=self.connected[id1].union(self.connected[id2])
                        for other in combined_connected:
                            self.connected[other]=combined_connected
            #compare to parent cells
            for cell_body2, cell_geom2 in self.cellList[:newCellStart]:
                data2=cell_body2.getData()
                state2=data2["type"]
                if self.InteractionsNetwork[state1][state2]>0:
                    #pdb.set_trace()
                    id2=data2["id"]
                    if id2 not in parentStates:
                        parentStates[id2]=data2["type"]
                    ##################################################################################################################
                    ##### If the parent cells were joined                                                               ##############
                    ##################################################################################################################
                    if self.JointsNetwork[parentIDs[id1]][id2]>0:
                        pos1=cell_body.getPosition()
                        pos2=cell_body2.getPosition()
                        distance=((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)**.5
                        if distance<=2.2*radius:
                            j = OdeBallJoint(self.odeWorld)
                            j.attach(cell_body, cell_body2)
                            if id1 not in self.jointLookup:
                                self.jointLookup[id1]=dict()
                            if id2 not in self.jointLookup:
                                self.jointLookup[id2]=dict()
                            self.jointLookup[id1][id2]=(j, self.InteractionsNetwork[state1][state2])
                            self.jointLookup[id2][id1]=(j, self.InteractionsNetwork[state1][state2])
                            if id1<id2:
                                self.noDetach[(id1,id2)]=0
                            else:
                                self.noDetach[(id2,id1)]=0
                            #j.setAnchor( (1,2,0) )
                            #print("Adding joint",(id1,id2))
                            self.JointsNetwork[id1][id2] = True
                            self.JointsNetwork[id2][id1] = True
                            #print("Adding joint")
                            if id2 not in self.connected[id1]:
                                combined_connected=self.connected[id1].union(self.connected[id2])
                                for other in combined_connected:
                                    self.connected[other]=combined_connected
            #compare to daughter cells but do same comparison of parents, make sure not to do comparison twice.
            for cell_body2, cell_geom2 in self.cellList[currentIndex+1:]:
                data2=cell_body2.getData()
                state2=data2["type"]
                if self.InteractionsNetwork[state1][state2]>0:
                    #pdb.set_trace()
                    id2=data2["id"]
                    if id2 not in parentStates:
                        parentStates[id2]=data2["type"]
                    ##################################################################################################################
                    ##### If the parent cells were joined                                                               ##############
                    ##################################################################################################################
                    if self.JointsNetwork[parentIDs[id1]][parentIDs[id2]]>0:
                        pos1=cell_body.getPosition()
                        pos2=cell_body2.getPosition()
                        distance=((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)**.5
                        if distance<=2.2*radius:
                            j = OdeBallJoint(self.odeWorld)
                            j.attach(cell_body, cell_body2)
                            if id1 not in self.jointLookup:
                                self.jointLookup[id1]=dict()
                            if id2 not in self.jointLookup:
                                self.jointLookup[id2]=dict()
                            self.jointLookup[id1][id2]=(j, self.InteractionsNetwork[state1][state2])
                            self.jointLookup[id2][id1]=(j, self.InteractionsNetwork[state1][state2])
                            if id1<id2:
                                self.noDetach[(id1,id2)]=0
                            else:
                                self.noDetach[(id2,id1)]=0
                            #j.setAnchor( (1,2,0) )
                            #print("Adding joint",(id1,id2))
                            self.JointsNetwork[id1][id2] = True
                            self.JointsNetwork[id2][id1] = True
                            #print("Adding joint")
                            if id2 not in self.connected[id1]:
                                combined_connected=self.connected[id1].union(self.connected[id2])
                                for other in combined_connected:
                                    self.connected[other]=combined_connected
            currentIndex+=1
        ##################################################################################################################
        ##### Checks if cells were connected and should still be connected ###############################################
        ##################################################################################################################
        for cell_body, cell_geom in self.cellList[:newCellStart]:
            data1=cell_body.getData()
            id1=data1["id"]
            state1=data1["type"]
            if id1 in self.jointLookup:
                detached=[]
                for id2 in self.jointLookup[id1]:
                    if id2>id1:
                        body2=self.cellList[id2][0]
                        data2=body2.getData()
                        state2=data2["type"]
                        if self.InteractionsNetwork[state1][state2]<=0:
                            joint=self.jointLookup[id1][id2][0]
                            joint.detach()
                            self.JointsNetwork[id1][id2] = False
                            self.JointsNetwork[id2][id1] = False
                            detached.append(id2)
                        #updates bond strength if it has changed
                        elif self.InteractionsNetwork[state1][state2]!=self.jointLookup[id1][id2][1]:
                            joint=self.jointLookup[id1][id2][0]
                            self.jointLookup[id1][id2]=(joint,self.InteractionsNetwork[state1][state2])
                            self.jointLookup[id2][id1]=(joint,self.InteractionsNetwork[state1][state2])
                for id2 in detached:
                    self.jointLookup[id1].pop(id2)
                    self.jointLookup[id2].pop(id1)
                    self.fixConnected(id1, id2)
                
    

                                
    def updateJoints3(self, cellList, newCellStart, parentStates, parentIDs):
        radius=6
        oldInteractionsNetwork=self.InteractionsNetwork
        self.InteractionsNetwork=self.FSM.add_interaction_energies_fast(cellList)
        for cell_body,cell_geom in cellList:
            data1=cell_body.getData()
            id1=data1["id"]
            state1=data1["type"]
            if id1<newCellStart:
                if id1 in self.jointLookup:
                    for id2 in self.jointLookup[id1]:
                        joint=self.jointLookup[id1][id2]
                        body2=joint[5]
                        data2=body2.getData()
                        id2=data2["id"]
                        state2=data2["type"]
                        if not self.InteractionsNetwork[state1][state2]:
                            joint[0].detach()
                            joint[1].detach()
                            joint[2].detach()
                            joint[0].destroy()
                            joint[1].destroy()
                            joint[2].destroy()
                            joint[3].destroy()
                            joint[4].destroy()
                            self.JointsNetwork[id1][id2] = False
                            self.JointsNetwork[id2][id1] = False
                            self.jointLookup[id1].pop(id2)
            else:
                for cell_body2, cell_geom2 in cellList[:newCellStart]:
                    data2=cell_body2.getData()
                    state2=data2["type"]
                    if self.InteractionsNetwork[state1][state2]:
                        #pdb.set_trace()
                        id2=data2["id"]
                        if oldInteractionsNetwork[parentStates[parentIDs[id1]]][parentStates[id2]]:
                            pos1=cell_body.getPosition()
                            pos2=cell_body2.getPosition()
                            distance=((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)**.5
                            if distance<=2.2*radius:
                                self.create_3Part_Joint(cell_body, cell_body2, id1, id2)
                                #j.setAnchor( (1,2,0) )
                                #print("Adding joint",(id1,id2))
                                self.JointsNetwork[id1][id2] = True
                                self.JointsNetwork[id2][id1] = True
                                #print("Adding joint")
            
    def mitosis(self):
        """
        Function: Handles functionality realted to mitosis/cell division and state changes.
    
        Steps:
        Increment cell division counter. Perform cell division if counter over threshold.
        Perform position corrections (destroy joints, remove force/velocity, mass spring adjustment, rebuild joints)
    
        Temp shortcut to cell data
        myBody.setData({"id":str(cell_id),"type":0,"mitosis_trigger":20,"mitosis_counter":0})
        """
    
        # perorm cell divisions. cytokenesis
        mitosis_triggered = False
        firstNewCell=True
        firstNewCellID=0
        parentIDs=dict()
        parentStates=dict()
        for idx in range(len(self.cellList)):
            cell_body,cell_geom = self.cellList[idx] #cellList is a list of tuples.
            
            #rendering
            #cell_geom.setPosQuat(cell_body.getPosition(), Quat(cell_body.getQuaternion()))
            x,y,z = cell_body.getPosition()
            u,v,w = cell_body.getLinearVel()
            cell_data = cell_body.getData()
            cell_data["mitosis_counter"] += self.dt
    
            if cell_data["mitosis_trigger"] >= 0 and cell_data["mitosis_counter"] >= cell_data["mitosis_trigger"]:
                #print("Adding new cell")
                mitosis_triggered = True
    
                parent_cell_body = cell_body
                parent_cell_data = cell_data
                parent_state_id = parent_cell_data["type"]
                
                # Handle Asymmetric Cell Division.
                state1, state2 = self.FSM.getCellTransitions(cell_state_id=parent_state_id)
    
                # Update parent(reset clocks, change state_id,check cell  cycle arrest)
                new_parent_data = self.FSM.updateAttributes(parent_cell_data,new_state_id=state1)
                parent_cell_body.setData(new_parent_data)
    
                # Create and update child cell.
                # place new cell with new body and geom. give parent attributes
                if firstNewCell:
                    firstNewCell=False
                    firstNewCellID=len(self.cellList)
                child_id = len(self.cellList)
                parentIDs[child_id]=parent_cell_data["id"]
                parentStates[parent_cell_data["id"]]=parent_state_id
                parent_position = np.array([x,y,z])
                half_radius = 3.0
                child_position = self.RandomPointOnSphere()*half_radius + parent_position
                child_cell_obj = self.create_cell(self.odeWorld,
                                        position=child_position,
                                        init_force=[0,0,0], #no initial force
                                        cell_id=child_id,
                                        cell_type=state2)
                # Update child attributes
                child_body = child_cell_obj[0]
                child_data = child_body.getData()
                new_child_data = self.FSM.updateAttributes(child_data,new_state_id=state2)
                child_body.setData(new_child_data)
    
                self.cellList.append(child_cell_obj)
    
                
                #process state change. updates mitosis trigger, type,
        ##############################################
        # Reposition cells after mitosis.
        ############################################
        if mitosis_triggered == True:
            # remove joints
            #self.remove_joints(self.cellList)
            #self.JointsNetwork = self.createJointNetwork()
            #TD change, update joints, instead of always breaking
            self.updateJoints(firstNewCellID, parentStates, parentIDs)
            #adjust positions
            #cellList = OpitmizeCellPositions(cellList, self.JointsNetwork, self.optimizationTime)
            #SimulationAnalysis.report_distances(cellList)
    
            # recreate joints.
            #This will be handled in next update step.
    
            #create allowed interactions network with current cell population
            
            #TD not needed becasue now done in updateJoints:
            #self.InteractionsNetwork = self.FSM.add_interaction_energies_fast(cellList)
    
    def remove_joints(self):
        for cell_body,cell_geom in self.cellList:
            for joint in cell_body.getJoints():
                joint.detach()

if __name__=="__main__":
    if 0==1:
        simulator=Run_Simulation()
        simulator.begin("Lim_Green_Red","examples/stateMachine_Lim_Green_Red_002_max_simple.csv")
    if 1==0:
        simulator=Run_Simulation()
        simulator.begin("Letter_D","examples/stateMachine8count_7tet_lowerD_numb_early_002_max_simple.csv")
    if 0==1:
        simulator=Run_Simulation()
        simulator.begin("ape1_ape2","examples/stateMachine_ape1_ape2_test.csv")
    if 0==1:
        simulator=Run_Simulation()
        simulator.begin("lim_green_red_new","examples/2018_01_18_lim_green_red_state_machine_simple.csv")
    if 0==1:
        simulator=Run_Simulation()
        simulator.begin("example_rod_max", "examples/test_rod_max_state_machine_002_simple.csv")
    if 0==1:
        simulator=Run_Simulation()
        simulator.begin("test_rod_max", "examples/test_rod_max_state_machine_002_simple.csv")
    if 0==1:
        block_reference=dict()
        block_reference["A"]=(0,False)
        block_reference["D"]=(1,False)
        block_reference["G"]=(2,False)
        simulator=Run_Simulation()
        simulator.begin(block_reference, "3tet_max", "examples/3tet_state_machine_002_simple.csv", 5)
    if 0==1:
        simulator=Run_Simulation()
        simulator.begin("ape1_6_large", "examples/stateMachine_ape1_6_v3_test.csv")
    if 1==1:
        block_reference=dict()
        block_reference["A"]=(0,False)
        block_reference["C"]=(1,False)
        block_reference["D"]=(2,False)
        simulator=Run_Simulation()
        simulator.begin(block_reference, "2019_06_04_3tet_high", "examples/2019_06_04_3tets_register_state_machine_simple.csv", 7)