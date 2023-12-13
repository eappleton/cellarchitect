"""
Author: Demarcus Briers, Boston University
    BU Robotics Lab, Hybrid and Networked Systems Group
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
import os,sys
import errno
from datetime import datetime
import time
import gzip
import pdb
import numpy as np
import shutil
from collections import defaultdict
from xml.dom import minidom
import random
import pdb


class FiniteStateMachine():

    def __init__(self,file, bindingStrengths):
        """ Initialization Function
        file - path to finite state machine file.

        """
        self.bindingStrengths=bindingStrengths
        state_machine = file
        
        #1. Create sim dir and copy simulation template to temp directory
        #src_path = "C:\\src\\CC3D\\simulations\\stateMachineSimulator"
        #out_path = "C:\\src\\CC3D\\simulations\\stateMachineOutput"
        #sim_folder = generate_sim_folder(src_path,out_path,state_machine)
        
        # 1. Read FSM
        self.states,self.generic_binders,self.paired_binders,self.transitions, self.expressedIgnored = self.read_states_and_transitions(state_machine)
        #self.states = states
        #self.transitions = transitions

        #2. Create jointnetwork. record of connected cells.
        #self.jointsnetwork = createJointNetwork(max_cells=200)
        
        #3. Create interactions/adhesion network.
        #self.interactions_network = self.add_interaction_energies()

        

    def generate_sim_folder(self,state_machine,project_label,output_folder="simulations_1.0"):
        """ Generate sim folder using current time stamp """
        
        #make simulation directory
        time_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        sim_folder = os.path.join(
            output_folder, 
            project_label+"_"+time_stamp)
        
        # copy template files to sim directory.
        try:
            os.makedirs(sim_folder)
            #shutil.copytree(src_folder,sim_folder)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise  # This was not a "directory exist" error..
        
        #copy state machine so cell division mechanism can load from the file
        state_machine_copy = os.path.join(sim_folder,"state_machine.csv")
        #shutil.copy2(state_machine,state_machine_copy)
        
        print("Simulation template copied to %s" % sim_folder)
        return sim_folder

    def read_states_and_transitions(self,filename):
        """
        Function: parse a file that contains both states and transitions for
        simulations of 3D cell geometries. File can be plain CSV or gzipped CSV.
        Note: the file must end with .csv or .gz
        
        First N lines describe the cell state: id protein1;protein2;proteinN
        follow by a blank line,
        followd by M lines describing transition probabilities after mitosis:
            1 1 .05
            1 2
        """
        states = {}
        generic_binders = {}
        paired_binders = {}
        expressedIgnored= dict()
        transitions = {"transitions": defaultdict(lambda:[]),
                       "probabilities":defaultdict(lambda:[])}
        read_mode = "states" #states or transitions
        
        
        #--------------------#
        #       Read file    #
        #--------------------#
        with open(filename, "r") as f:
            header = f.readline()
            for line in f:
                line = line.strip()
                
                # check with read mode to be in
                if line == "":
                    read_mode = "transitions"
                    continue
                
                #process line
                if read_mode == "states":
                    states,generic_binders,paired_binders, expressedIgnored = self.read_state(states,
                                                                       generic_binders,
                                                                       paired_binders,
                                                                       line, expressedIgnored)
                elif read_mode == "transitions":
                    transitions = self.read_transition(transitions,line)
                else:
                    assert 0==1,"Unkown read mode %s " % (str(read_mode))
        #Assign random strengths to interactions with proteins not in the lookup file.
        all_generic=set()
        all_paired=set()
        for state in generic_binders:
            all_generic=all_generic.union(generic_binders[state])
            all_paired=all_paired.union(paired_binders[state])
            
        unspecified_generic=[]
        for protein in all_generic:
            if protein not in self.bindingStrengths:
                self.bindingStrengths[protein]=dict()
                unspecified_generic.append(protein)
        unspecified_paired=[]
        for protein in all_paired:
            if protein not in self.bindingStrengths:
                self.bindingStrengths[protein]=dict()
                unspecified_paired.append(protein)
        alreadyDone=set()
        for protein in unspecified_generic:
            for protein2 in self.bindingStrengths:
                if protein2 not in alreadyDone:
                    if protein!=protein2:
                        strength=0.3*random.random()
                        self.bindingStrengths[protein][protein2]=strength
                        self.bindingStrengths[protein2][protein]=strength
                    else:
                        self.bindingStrengths[protein][protein]=0.2*random.random()+.8
        for protein in unspecified_paired:
            for protein2 in self.bindingStrengths:
                if protein2 not in alreadyDone:
                    if protein2 in all_paired and protein[:-1]==protein2[:-1]:
                        strength=0.2*random.random()+.8
                        self.bindingStrengths[protein][protein2]=strength
                        self.bindingStrengths[protein2][protein]=strength
                    else:
                        strength=0.3*random.random()
                        self.bindingStrengths[protein][protein2]=strength
                        self.bindingStrengths[protein2][protein]=strength
        # update transitions to sum to 1.0
        for cell_id in transitions["transitions"].keys():
            probabilities_sum = np.sum(transitions["probabilities"][cell_id])
            
            remaining_other_prob = 1.0 - probabilities_sum
            if remaining_other_prob > 0:
                transitions["probabilities"][cell_id].append(remaining_other_prob)
                transitions["transitions"][cell_id].append("-1,-1")

        return states,generic_binders,paired_binders,transitions, expressedIgnored

    def read_state(self,states,generic_binders,paired_binders,line, expressed_ignored):
        """ Parse the mappings of state to surface proteins.

            Returns:
            states - mapping from state_id to expressed proteins (includes CellCycleArrest)
            generic_binders - mapping from state_id to generic surface binding proteins.
            pared_binders - mapping from state_id to paired surface binding proteins
        """
        ignored_proteins = {"iRFP670", "eGFP", "mScarlet", "tdTomato", "EBFP2", "EYFP", "mBanana", "CellDeath", "mPlum"}
        state_id, proteins_list = line.strip(";").split(",")  #omit trailing semicolon
        #proteins_list_raw = []
        
        # parse allowed proteins
        state_id=int(state_id)
        if proteins_list == "":
            proteins_list = set()
            expressed_ignored[state_id]=set()
        else:
            proteins_list = proteins_list.split(";")
            proteins_list = set(proteins_list)
            expressed_ignored[state_id]=proteins_list.intersection(ignored_proteins)
            proteins_list = proteins_list.difference(ignored_proteins)
        
        states[state_id] = proteins_list
        generic_binders[state_id]  = \
            set([x for x in proteins_list if not x.endswith('1') and not x.endswith('2') and x != "CycleArrest"])
        
        paired_binders[state_id]  = \
            set([x for x in proteins_list if x.endswith('1') or x.endswith('2')])
        #print(state_id)
        #print(proteins_list)
        
        return states,generic_binders,paired_binders, expressed_ignored

    def read_transition(self,transitions,line):
        """ Parse the mappings of state transition to probability
        Format: parent_state,child1_state,child2_state,transition_probability
        """
        #print(line.split(","))
        parent,child1,child2,probability = line.strip("").split(",")[:4]
        parent=int(parent)
        child_transitions = (int(child1),int(child2))
        transitions["transitions"][parent].append(child_transitions)  #paired chld cell states 
        transitions["probabilities"][parent].append(float(probability))  # transition probabilites
        
        return transitions

    def add_interaction_energies_fast(self,cellList):
        """
        Functon: Create graph that maps if pairs of cell_staes will adhere (having interacting surface proteins)
        This function calculates possible interactions between states currently in the simulation.
        This is much more effiecient that add_interaction_energies()

        Inputs:
        -----------
        cellList - tuple of current cell objects in the simulations. We extract the cell state ids from this to compute runtime
        interactions. This function scales better for relatively small geometric shapes <= 10,000 cell. 
        """
        state_ids = set()
        for cell_obj,cell_gom in cellList:
            cell_data = cell_obj.getData()
            cell_state_id = cell_data["type"]
            state_ids.add(cell_state_id)
        
        #num_states = len(state_ids)
        #N_interactions = (num_states**2 + num_states)/2  #Nth Triangular Number

        
        # Create interaction Network.
        I = defaultdict(lambda:{}) # Interaction Network.
        
        
        ######################################
        #Compute static cell-cell interactions
        ######################################
        #print_freq = np.round(0.1*num_states) #freqeucny of progress bar update
        #iteration = 0
        #print("Computing interactions for %i cell states and %i interactions" % (num_states,int(N_interactions)))
        #start = time.time()
        # interactions between different states
        for i in state_ids:
            for j in state_ids:
                if i>=j:
                    self.update_interaction_network(i,j,I)
                    #iteration+=1 #update the counter
                    #if iteration % print_freq == 0:
                    #    self.printProgressBar(iteration,int(N_interactions))
        
                   
        #end = time.time()
        #print("\nComputed %i Interactions in %s seconds" % (iteration,end-start))
        #print(I)
        # Return Interaction Netorks (vertices=state_id,edges=)
        #sys.exit()
        #print(self.bindingStrengths)
        #print(I)
        return I
    
    #deprecated
    def add_interaction_energies(self):
        """
        Functon: Create graph that maps if pairs of cell_staes will adhere (having interacting surface proteins)
        Note: this function exhaustively calculates all interactions between all possible states. A much smarte
        apporach is to calculate possible interactions between states currently in the simulation.
        """
        state_ids = sorted(self.states.keys(),key=int)
        num_states = len(self.states.keys())
        N_interactions = (num_states**2 + num_states)/2  #Nth Triangular Number

        #Check formatting of the file.
        state_ints = [int(x) for x in state_ids]
        assert max(state_ints) == num_states-1, "State machine state IDs are not zero indexed. This is an assumption for fast iterations."
        
        # Create interaction Network.
        I = defaultdict(lambda:{}) # Interaction Network.
        I = np.zeros((num_states,num_states),dtype=np.bool_)
        print_freq = np.round(0.1*num_states)
        iteration = 0
        #print("Computing interactions for %i cell states and %i interactions" % (num_states,int(N_interactions)))
        
        ######################################
        #Compute static cell-cell interactions
        ######################################
        start = time.time()
        # interactions between different states
        for i in range(0,num_states):
            for j in range(i+1,num_states):
                self.update_interaction_network(i,j,I)
                iteration+=1 #update the counter
                if iteration % print_freq == 0:
                    self.printProgressBar(iteration,int(N_interactions))
        
        # self interactions
        for i in range(0,num_states):
            j = i
            self.update_interaction_network(i,j,I)
            iteration+=1 #update the counter
            if iteration % print_freq == 0:
                    self.printProgressBar(iteration,int(N_interactions))
                   
        end = time.time()
        #print("\nComputed %i Interactions in %s seconds" % (iteration,end-start))
        #print(I["10"])
        # Return Interaction Netorks (vertices=state_id,edges=)
        #sys.exit()
        #print(I)
        return I

    def update_interaction_network(self,i,j,I):
        """
        i - index 1 (int)
        j - index 2 (int)
        I - interaction network, can be represent as numpy matrix, nested defaultdict, or netorkX network.
        """
        #if either cell is expressing no proteins, then we set binding affinity to -1 to indicate no attraction.
        num_i=len(self.generic_binders[i])+len(self.paired_binders[i])
        num_j=len(self.generic_binders[j])+len(self.paired_binders[j])
        if num_i==0 or num_j==0:
            I[i][j] = -1.0
            I[j][i] = -1.0
        else:
            strength=0.0
            #do all pairwise evaluations between the proteins of state i and state j.
            for protein in self.generic_binders[i]:
                for protein2 in self.generic_binders[j]:
                    strength+=self.bindingStrengths[protein][protein2]
                for protein2 in self.paired_binders[j]:
                    strength+=self.bindingStrengths[protein][protein2]
            for protein in self.paired_binders[i]:
                for protein2 in self.generic_binders[j]:
                    strength+=self.bindingStrengths[protein][protein2]
                for protein2 in self.paired_binders[j]:
                    strength+=self.bindingStrengths[protein][protein2]
            strength=strength/num_i/num_j #temporarily turned off averaging for BU sims.
            I[i][j] = strength
            I[j][i] = strength
        
    def getCellTransitions(self,cell_state_id):
        """
        Determine stochastic cell state transitions using the Finite State Machine. 
        """
        parent_id = cell_state_id

        states = self.transitions["transitions"][parent_id]
        probs = self.transitions["probabilities"][parent_id]

        # get transition states
        state1,state2 = random.choices(states, weights=probs)[0]

        #print("Division triggered.\nChild cell states %s" % child_cell_states)

        child_states = [state1,state2]
        np.random.shuffle(child_states)

        return child_states

    def updateAttributes(self,cell_data,new_state_id):
        """
        Update cell data after mitosis:
        1. reset the cell division clock
        2. update cell state id
        3. check if cellcyclearrest occurs

        Returns:
        cell_data dictionary
        """
        cell_data["mitosis_counter"] = 0
        cell_data["type"] = new_state_id

        #Determine cell cycle arrest
        if "CycleArrest" in self.states[new_state_id]:
            cell_data["mitosis_trigger"] = -1 #indicated no cell division. See run_simulation.mitosis() also

        return cell_data


    def printProgressBar (self,iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = ' '):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
        # Print New Line on Complete
        if iteration == total: 
            print()
