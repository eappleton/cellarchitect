import itertools
import csv
import math
import pdb
import StateMachineDrawings as draw
import time
from multiprocessing import Lock, Pipe, Process, Queue, Array, Value
import multiprocessing
import queue
import inspect
import os, psutil

class StateMachine():
    
    def __init__(self, outformat, outFile, counterFile, registerFile=None, recLevel="High", unitTest=False):
        '''
        Constructor.
        
        Authors:
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 9/11/2018
        
        Arguments:
        outformat -- list of formats that should be outputted
        outFile -- path for the "Detailed" output file.
        counterFile -- path for the csv containing the counter circuit.
        registerFile -- path for the csv containing the register circuit.
        recLevel -- "Low", "Mid", "High", or "Max". Indicates which recombinase file to use and which
            numb probability values to use.
        
        Efficiency:
        O(n^2) where n is the total number of circuit components.
        '''
        if not unitTest:
            self.outFile=outFile
            
            #splitNumb, bothNumb
            if recLevel!="Max":
                self.numbValues={0: .85, 1: .15}
            else:
                self.numbValues={0: 1, 1: 0}
            #chance that if 1 recombinase is correctly segregated that another will be.
            self.numbCorrelation=[.75,.25]
            self.nextTimeStamp=time.time()+60
            self.outformat=outformat
            self.pdf=outFile[:-3]+"pdf"
            
            #initialize data structures
            self.arrestID=set()
            self.textToInt=dict()
            self.intToText={0:""}
            self.allCassettes=dict()
            self.allRDFs=dict()
            self.allSiteSets=dict()
            self.allPromoters=set()
            self.allTerminators=set()
            self.terminatorIDs=set()
            self.numbIDs=dict()
            self.siteToRec=dict()
            self.recToRDF=dict()
            #will hold all edges between configurations.
            self.allEdges=[]
            #the index of the initial state. Will be updated as necessary
            self.startIndex=0
            #will hold all the found configurations
            self.allConfigs=dict()
            
            #reads counter from file
            reader = csv.reader(open(counterFile))
            self.counter = [x[0] for x in list(reader)]
            #represents circuit as list of ints
            self.identifyComponents(self.counter)
            self.intCounter = self.encode(self.counter)
            
            #reads register from file
            if registerFile:
                reader = csv.reader(open(registerFile))
                self.register= [x[0] for x in list(reader)]
                #represents circuit as list of ints
                self.identifyComponents(self.register)
                self.intRegister = self.encode(self.register)
                #adds the initial configuration
                self.allConfigs[(tuple(self.intCounter),tuple(self.intRegister))]=0
            #deals with cases where there is no register
            else:
                self.register=None
                self.intRegister=None
                #adds the initial configuration
                self.allConfigs[(tuple(self.intCounter),)]=0
                            
            #maps the code for each site to the code for its recombinase
            self.mapSites()
            
            #checks to see if the constructor has been called from within this module or from Unite.py, because this changes the path
            #needed to find the recombinase efficiencies file.
            pyFile=inspect.stack()[-1].filename
            
            #reads recombinase efficiencies from file
            if pyFile=="StateMachine.py":
                if recLevel=="High":
                    reader=csv.reader(open("rec_efficiencies/Recombinase_Efficiencies_High.csv"))
                elif recLevel=="Mid":
                    reader=csv.reader(open("rec_efficiencies/Recombinase_Efficiencies_Mid.csv"))
                elif recLevel=="Low":
                    reader=csv.reader(open("rec_efficiencies/Recombinase_Efficiencies_Low.csv"))
                else:
                    reader=csv.reader(open("rec_efficiencies/Recombinase_Efficiencies_Max.csv"))
            #if the constructor has been called from Unite.py
            else:
                if recLevel=="High":
                    reader=csv.reader(open("circuit_design/rec_efficiencies/Recombinase_Efficiencies_High.csv"))
                elif recLevel=="Mid":
                    reader=csv.reader(open("circuit_design/rec_efficiencies/Recombinase_Efficiencies_Mid.csv"))
                elif recLevel=="Low":
                    reader=csv.reader(open("circuit_design/rec_efficiencies/Recombinase_Efficiencies_Low.csv"))
                else:
                    reader=csv.reader(open("circuit_design/rec_efficiencies/Recombinase_Efficiencies_Max.csv"))
            
            #reads in rec efficiencies
            self.recEfficiency=dict()
            firstRecProb=[]
            for x in list(reader)[1:]:
                fwd=x[0]+"_cassette_F"
                if not firstRecProb:
                    firstRecProb=[float(y) for y in x[1:]]
                if fwd in self.textToInt:
                    probs=[float(y) for y in x[1:]]
                    self.recEfficiency[self.textToInt[fwd]]=probs
                    self.recEfficiency[self.textToInt[fwd]-2]=probs
            #assigns rec efficiency to generic recombinase names.
            for x in self.allCassettes:
                if "Rec" in x and "Numb" not in x:
                    self.recEfficiency[self.allCassettes[x]]=firstRecProb
                    self.recEfficiency[self.allCassettes[x]+2]=firstRecProb
                    print(x+" assigned first probability.")
            
            #set of recombinases (necessary for differentiating recombinases from other expressed genes for simple outputs)
            self.recombinases={"int2", "int3", "int4", "int5", "int7", "int8", "int9", "int10", "int11", "int12", "int13",
                               "bxb1", "phic31", "tp901", "phirv1", "phibt1", "phijoe", "spra"}
            for i in range(len(self.counter)+len(self.register)):
                self.recombinases.add("rec"+str(i))
    
    
    def mapSites(self):
        '''
        Maps the code for each recombinase site to the code for its recombinase. Also maps each recombinase
        to its RDF.
        
        Authors:
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 9/11/2018
        
        Efficiency:
        O(n^2*m) where n=# of distinct circuit components and m is the longest recombinase name. Because the length of the
        recombinase names are generally pretty consistent, this is essentially O(n^2).
        '''
        for name in self.allSiteSets:
            for name2 in self.allCassettes:
                if name2==name[:name.index("_")]:
                    self.siteToRec[self.allSiteSets[name]]=self.allCassettes[name2]
                    break
        for name in self.allCassettes:
            for name2 in self.allRDFs:
                if name==name2:
                    self.recToRDF[self.allCassettes[name]]=self.allRDFs[name2]
                    break
                    
    
    def encode(self, sequence):
        '''
        Translates string list version of a circuit to int list.
        
        Authors:
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 9/11/2018
        
        Arguments:
        sequence -- list of strings that represents a circuit
        
        Efficiency:
        O(n) where n is the length of the list.
        '''
        encoded=[]
        for item in sequence:
            encoded.append(self.textToInt[item])
        return encoded
    
    def decode(self, sequence):
        '''
        Translates int list version of a circuit to string list.
        
        Authors:
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 9/11/2018
        
        Arguments:
        sequence -- list of ints that represents a circuit
        
        Efficiency:
        O(n) where n is the length of the list.
        '''
        decoded=[]
        for item in sequence:
            decoded.append(self.intToText[item])
        return decoded
    
    def identifyComponents(self, sequence):
        '''
        Determines a int representation of each compoment in a genetic circuit.
        
        Authors:
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 9/11/2018
        
        Arguments:
        sequence -- list of strings that represents a circuit
        
        Efficiency:
        O(n*m) where n is the total number of components and m is the longest component name. Since m is usually fairly consistent,
        this is essentially O(n)
        
        Maps each promoter and terminator to negative even numbers, the sites for each recombinase to positive even numbers,
        the genes that are not RDFs to positive odd numbers, and the RDFs to negative odd numbers. Fwd version of component is
        always code+2 relative to negative version.
        '''
        #pdb.set_trace()
        #finds the highest negative even number that hasn't been used yet
        pt=-2
        while pt in self.intToText:
            pt+=-4
        #finds the lowest positive even number that hasn't been used yet
        s=16
        while s in self.intToText:
            s+=16
        #finds the lowest positive odd number that hasn't been used yet
        c=1
        while c in self.intToText:
            c+=4
        #finds the highest negative odd number that hasn't been used yet
        r=-1
        while r in self.intToText:
            r+=-4
        for item in sequence:
            first=item.index("_")
            if "terminator" in item:
                name=item[:first]
                if name not in self.allTerminators:
                    #fwd version
                    self.textToInt[name+"_terminator_F"]=pt
                    self.intToText[pt]=name+"_terminator_F"
                    #rev version
                    self.textToInt[name+"_terminator_R"]=pt-2
                    self.intToText[pt-2]=name+"_terminator_R"
                    self.allTerminators.add(name)
                    self.terminatorIDs.add(pt)
                    self.terminatorIDs.add(pt-2)
                    pt+=-4
            elif "promoter" in item:
                name=item[:first]
                if name not in self.allPromoters:
                    #fwd version
                    self.textToInt[name+"_promoter_F"]=pt
                    self.intToText[pt]=name+"_promoter_F"
                    #rev version
                    self.textToInt[name+"_promoter_R"]=pt-2
                    self.intToText[pt-2]=name+"_promoter_R"
                    self.allPromoters.add(name)
                    pt+=-4
            elif "cassette" in item:
                name=item[:first]
                if name not in self.allCassettes:
                    #fwd version
                    self.textToInt[name+"_cassette_F"]=c+2
                    self.intToText[c+2]=name+"_cassette_F"
                    #rev version
                    self.textToInt[name+"_cassette_R"]=c
                    self.intToText[c]=name+"_cassette_R"
                    self.allCassettes[name]=c
                    if name=="CycleArrest":
                        self.arrestID={c,c+2}
                    #if the gene is mapped to numb, maps the numb fused version to the normal version
                    elif "Numb" in item:
                        rec=name[:name.index("-")]
                        if rec+"_cassette_R" in self.textToInt:
                            self.numbIDs[c]=self.textToInt[rec+"_cassette_R"]
                            self.numbIDs[c+2]=self.textToInt[rec+"_cassette_R"]
                        else:
                            self.numbIDs[c]=c+4
                            self.numbIDs[c+2]=c+6
                            self.textToInt[rec+"_cassette_F"]=c+6
                            self.intToText[c+6]=rec+"_cassette_F"
                            self.textToInt[rec+"_cassette_R"]=c+4
                            self.intToText[c+4]=rec+"_cassette_R"
                            self.allCassettes[rec]=c+4
                            c+=4
                    c+=4
            elif "RDF" in item:
                name=item[:first]
                if name not in self.allRDFs:
                    #fwd version
                    self.textToInt[name+"_RDF_F"]=r
                    self.intToText[r]=name+"_RDF_F"
                    #rev version
                    self.textToInt[name+"_RDF_R"]=r-2
                    self.intToText[r-2]=name+"_RDF_R"
                    self.allRDFs[name]=r-2
                    r+=-4
            #all sites are mapped to a range of 16.
            elif "site" in item:
                name=item[:item.index("_",first+1)]
                if name not in self.allSiteSets:
                    self.textToInt[name+"_attB_site_R"]=s
                    self.intToText[s]=name+"_attB_site_R"
                    self.textToInt[name+"_attB_site_F"]=s+2
                    self.intToText[s+2]=name+"_attB_site_F"
                    self.textToInt[name+"_attP_site_R"]=s+4
                    self.intToText[s+4]=name+"_attP_site_R"
                    self.textToInt[name+"_attP_site_F"]=s+6
                    self.intToText[s+6]=name+"_attP_site_F"
                    self.textToInt[name+"_attL_site_R"]=s+8
                    self.intToText[s+8]=name+"_attL_site_R"
                    self.textToInt[name+"_attL_site_F"]=s+10
                    self.intToText[s+10]=name+"_attL_site_F"
                    self.textToInt[name+"_attR_site_R"]=s+12
                    self.intToText[s+12]=name+"_attR_site_R"
                    self.textToInt[name+"_attR_site_F"]=s+14
                    self.intToText[s+14]=name+"_attR_site_F"
                    self.allSiteSets[name]=s
                    s+=16
            else:
                raise ValueError("Invalid Sequence Component: "+item)
    
    def getExpressedGenes(self, sequence, expressedCassettes, expressedRDFs):
        '''
        Identifies which genes are under promotion in a sequence.
        
        Authors:
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 9/11/2018
        
        Arguments:
        sequence -- list of ints that represents a circuit
        expressedCassettes -- set of expressed non-RDF genes.
        expressedRDFs -- set of expressed RDFs.
        
        Efficiency:
        O(n) where n is the length of the sequence
        
        Returns True if the cycle arrest gene is expressed, False otherwise.
        '''
        cycleArrest=False
        #checks forward strand
        prom=False
        for item in sequence:
            #if the component is a promoter, terminator, or RDF
            if item<0:
                #if the component is a forward promoter or terminator
                if item%4==2:
                    if item in self.terminatorIDs:
                        prom=False
                    else:
                        prom=True
                #if the component is under promotion and is a forward RDF.
                elif prom and item%4==3:
                    expressedRDFs.add(item-2)
            #if the component is under pomortion and is a non-RDF gene
            elif prom and item%4==3:
                expressedCassettes.add(item-2)
                if item in self.arrestID:
                    cycleArrest=True
        prom=False
        #checks reverse strand. same as above but for rev direction.
        for item in sequence[::-1]:
            if item<0:
                if item%4==0:
                    if item in self.terminatorIDs:
                        prom=False
                    else:
                        prom=True
                elif prom and item%4==1:
                    expressedRDFs.add(item)
            elif prom and item%4==1:
                expressedCassettes.add(item)
                if item in self.arrestID:
                    cycleArrest=True
        return cycleArrest
    
    def getActiveSites(self, expressedCassettes, expressedRDFs, sequence, windows, count, isCounter, recToWindows):
        '''
        Finds pairs of compatible recombinase sites whose recombinase is expressed.
        
        Authors: 
        Michaël Moret. Church Lab, Harvard Medical School. michael.moret@epfl.ch
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 1/8/2019
        
        Arguments:
        expressedCassettes -- set of int codes of the expressed genes
        expressedRDFs -- set of int codes of the expressed RDFs
        sequence -- list of int codes representing the elements in the DNA sequence.
        windows -- dictionary containing recombinase site pairs (windows) indexed by a pair id.
        count -- the next recombinase site pair id to use.
        isCounter -- True if the sequence is the counter, False if it is the register
        recToWindows -- Dictionary that contains list of rec site pair ids indexed by their recombinase's id
        
        Efficiency:
        O(n^2) where n is the length of the sequence
        
        Finds pairs of compatible recombinase sites whose recombinase is expressed, and puts them in windows, as well
        as adding them to the lookup recToWindows. Returns the next site pair id.
        '''
        allBPSites=dict()
        allLRSites=dict()
        #finds all the bp and lr sites in the sequence
        for i,x in enumerate(sequence):
            #if x is a site
            if x>0 and x%2==0:
                siteID=x-(x%16)
                if siteID in self.siteToRec and self.siteToRec[siteID] in expressedCassettes:
                    if x%16<7:
                        if self.siteToRec[siteID] in allBPSites:
                            if siteID in allBPSites[self.siteToRec[siteID]]:
                                allBPSites[self.siteToRec[siteID]][siteID].append((x,i))
                            else:
                                allBPSites[self.siteToRec[siteID]][siteID]=[(x,i)]
                        else:
                            allBPSites[self.siteToRec[siteID]]={siteID:[(x,i)]}
                    else:
                        if self.siteToRec[siteID] in allLRSites:
                            if siteID in allLRSites[self.siteToRec[siteID]]:
                                allLRSites[self.siteToRec[siteID]][siteID].append((x,i))
                            else:
                                allLRSites[self.siteToRec[siteID]][siteID]=[(x,i)]
                        else:
                            allLRSites[self.siteToRec[siteID]]={siteID:[(x,i)]}
                            
        for cassette in expressedCassettes:
            rdfActive=False
            if cassette in self.recToRDF and self.recToRDF[cassette] in expressedRDFs:
                rdfActive=True
            if cassette in allBPSites:
                bpSites = allBPSites[cassette]
            else:
                bpSites = []
            if cassette in allLRSites:
                lrSites = allLRSites[cassette]
            else:
                lrSites = []
            #Divies up the found sites into sets of orthogonal pairs. If at least one of each site in
            #the pair is found, the set is added to the return list.
            #does bp sites first
            for siteID in bpSites:
                for site1 in bpSites[siteID]:
                    if site1[0]%8<3:
                        for site2 in bpSites[siteID]:
                            if site2[0]%8>3:
                                #[start index, end index, is BP?, Recombinase ID, RDF active?, is counter?]
                                windows[count]=[min(site1[1], site2[1]), max(site1[1], site2[1]),True, cassette, rdfActive, isCounter]
                                recToWindows[cassette].append(count)
                                count+=1
            for siteID in lrSites:
                for site1 in lrSites[siteID]:
                    if site1[0]%8<3:
                        for site2 in lrSites[siteID]:
                            if site2[0]%8>3:
                                windows[count]=[min(site1[1], site2[1]), max(site1[1], site2[1]),False, cassette, rdfActive, isCounter]
                                recToWindows[cassette].append(count)
                                count+=1
        return count
             
    def branch(self, windows, conflicts, sequence, dependencies, combinations):
        '''
        Finds all combinations of windows that can occur simultaneously.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 9/11/2018
        
        Arguments:
        windows -- dictionary of windows, indexed by an ID for the window. Each window has form
            [start index, end index, is BP?, Recombinase ID, RDF active?, is counter?]
        conflicts -- A dictionary indexed by window IDs, containing a set of window IDs for each ID. The window
            IDs in the list are for windows that conflict with the indexing window.
        sequence -- The genetic sequence, represented as a list of strings.
        dependencies -- A list of tuples that keeps track of events that can only happen when another event doesn't happen,
            but not vice versa.
        combinations -- A dictionary of frozensets of combinations of window IDs
        
        Efficiency:
        Bad. Need to analyze closely.
    
        Finds all combinations of windows that can occur simultaneously. First identifies all the conflicts between windows,
        then creates buckets of conflicting combinations of windows, then finds all combinations involving at most 1 combination from each bucket.
        '''
        
        ##########################################################################
        ########## We identify conflicting windows ###############################   
        ##########################################################################  
        
        #pdb.set_trace()
        keys=list(windows.keys())
        conflicts2=dict()
        #For all pairs of windows. 
        for x, i in enumerate(keys):
            for j in keys[x+1:]:
                #share a site. mutually exclusive in all cases
                if windows[i][0]==windows[j][0] or windows[i][0]==windows[j][1] or windows[i][1]==windows[j][0] or windows[i][1]==windows[j][1]:
                    #updates the dictionary of conflicts.
                    if i in conflicts:
                        conflicts[i].add(j)
                    else:
                        conflicts[i]=set([j])
                    if j in conflicts:
                        conflicts[j].add(i)
                    else:
                        conflicts[j]=set([i])
                    if i in conflicts2:
                        conflicts2[i].add(j)
                    else:
                        conflicts2[i]=set([j])
                    if j in conflicts2:
                        conflicts2[j].add(i)
                    else:
                        conflicts2[j]=set([i])
                elif windows[i][0]<windows[j][0]<windows[i][1]:
                    #ABBA
                    if windows[i][0]<windows[j][1]<windows[i][1]:
                        #A excision. A excludes B
                        if (sequence[windows[i][0]]-sequence[windows[i][1]])%4==0:
                            if j in conflicts:
                                conflicts[j].add(i)
                            else:
                                conflicts[j]=set([i])
                            if i in conflicts2:
                                conflicts2[i].add(j)
                            else:
                                conflicts2[i]=set([j])
                            if j in conflicts2:
                                conflicts2[j].add(i)
                            else:
                                conflicts2[j]=set([i])
                            #adds to dependencies because B depends on A not happening.
                            dependencies.append((j,i))
                        #ELSE: A inversion. Both are fine
                    #ABAB. Always mutually exclusive
                    else:
                        if i in conflicts:
                            conflicts[i].add(j)
                        else:
                            conflicts[i]=set([j])
                        if j in conflicts:
                            conflicts[j].add(i)
                        else:
                            conflicts[j]=set([i])
                        if i in conflicts2:
                            conflicts2[i].add(j)
                        else:
                            conflicts2[i]=set([j])
                        if j in conflicts2:
                            conflicts2[j].add(i)
                        else:
                            conflicts2[j]=set([i])
                elif windows[j][0]<windows[i][0]<windows[j][1]:
                    #BAAB
                    if windows[j][0]<windows[i][1]<windows[j][1]:
                        #B exicision. B exludes A
                        if (sequence[windows[j][0]]-sequence[windows[j][1]])%4==0:
                            if i in conflicts:
                                conflicts[i].add(j)
                            else:
                                conflicts[i]=set([j])
                            if i in conflicts2:
                                conflicts2[i].add(j)
                            else:
                                conflicts2[i]=set([j])
                            if j in conflicts2:
                                conflicts2[j].add(i)
                            else:
                                conflicts2[j]=set([i])
                            #adds to dependencies because A depends on B not happening.
                            dependencies.append((i,j))
                        #ELSE: B inversion. Both are fine
                    #BABA. Always mutually exclusive
                    else:
                        if i in conflicts:
                            conflicts[i].add(j)
                        else:
                            conflicts[i]=set([j])
                        if j in conflicts:
                            conflicts[j].add(i)
                        else:
                            conflicts[j]=set([i])
                        if i in conflicts2:
                            conflicts2[i].add(j)
                        else:
                            conflicts2[i]=set([j])
                        if j in conflicts2:
                            conflicts2[j].add(i)
                        else:
                            conflicts2[j]=set([i])
                            
        ##########################################################################
        ########## We make buckets of conflicting windows ########################   
        ##########################################################################
        print("conflicts2",conflicts2)
        groups=[]
        done_windows=set()
        windowToBucket=dict()
        for key in keys:
            if key not in done_windows:
                windowToBucket[key]=len(groups)
                groups.append([key])
                done_windows.add(key)
                windowToBucket[key]=len(groups)-1
                self.getAllConflicts(groups[-1], key, done_windows, conflicts2, len(groups)-1, windowToBucket)
        buckets=[]
        print("groups", groups)
        for group in groups:
            buckets.append(self.createBucket(group, conflicts2))
        ##########################################################################
        ##### We find all combinations with only one group from each bucket  #####
        ##########################################################################  
        print("buckets",buckets)
        for bucket in buckets:
            combinations2=combinations.copy()
            for combo in combinations2:
                for item in bucket:
                    newCombo=set(combo)
                    for num in item:
                        newCombo.add(num)
                    combinations[frozenset(newCombo)]=0
        return buckets, windowToBucket
    
    def getAllConflicts(self, group, current, done_windows, conflicts2, bucket_id, windowToBucket):
        if current in conflicts2:
            for conflict in conflicts2[current]:
                if conflict not in done_windows:
                    done_windows.add(conflict)
                    group.append(conflict)
                    windowToBucket[conflict]=bucket_id
                    self.getAllConflicts(group, conflict, done_windows, conflicts2, bucket_id, windowToBucket)
                
    def createBucket(self, group, conflicts2):
        if len(group)>1:
            combinations=self.createBucket(group[:-1], conflicts2)
            combinations2=combinations.copy()
            for c in combinations2:
                conflicting=False
                for x in c:
                    if group[-1] in conflicts2[x]:
                        conflicting=True
                        break
                if not conflicting:
                    c2=c.copy()
                    c.append(group[-1])
                    combinations.append(c2)
            combinations.append([group[-1]])
            return combinations
        elif len(group)==1:
            return [[group[0]]]
        else:
            return []
                
    
    def event(self, window, combo, counter, register):
        '''
        Excises or Inverts a section of a circuit defined by window.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/23/2018
        
        Arguments:
        window -- [start index, end index, is BP?, Recombinase ID, RDF active?, is counter?]
        combo -- dictionary containing other windows
        counter -- the counter sequence, represented as a list of int codes
        register -- the register sequence, represented as a list of int codes. can be None
        
        Efficiency:
        O(n+m) where n is the length of the configuration and m is the size of combo.
    
        Excises or Inverts a section of a circuit defined by window. Updates the indices of the other windows in combo
        in the case of an inversion.
        '''
        #pdb.set_trace()
        if window[-1]:
            if counter[window[0]]%4==counter[window[1]]%4:
                self.excise(counter, window[0], window[1], window[2])
            else:
                self.invert(counter, window[0], window[1], window[2])
                self.updateWindows(window[0], window[1], combo, False)
        else:
            if register[window[0]]%4==register[window[1]]%4:
                self.excise(register, window[0], window[1], window[2])
            else:
                self.invert(register, window[0], window[1], window[2])
                self.updateWindows(window[0], window[1], combo, True)
    
    def updateWindows(self, bottom, top, combo, isRegister):
        '''
        Updates the indices of all recombinase windows within a region of a sequence to reflect an inversion.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/23/2018
        
        Arguments:
        bottom -- the lower index defining the range
        top -- the upper index defining the range
        combo -- dictionary containing the windows
        isRegister -- boolean of whether the inversion occured on the counter or the register.
        
        Efficiency:
        O(n) where n is the number of windows in combo
    
        Updates the indices of all recombinase windows within a region of a sequence to reflect an inversion.
        Checks if they are on the correct sequence (register or counter) and then changes their indices.
        '''
        #    pdb.set_trace()
        for i in combo:
            window=combo[i]
            if (isRegister and not window[-1]) or ((not isRegister) and window[-1]):
                temp=window.copy()
                if bottom<window[0]<top:
                    temp[0]=top-(window[1]-bottom)
                if bottom<window[1]<top:
                    temp[1]=top-(window[0]-bottom)
                combo[i]=temp

    def excise(self, sequence, first, last, isBP):
        '''
        Excises elements from a dna sequence by replacing them with 0.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/23/2018
        
        Arguments:
        sequence -- The dna sequence. A list of int codes
        first -- the lower index defining the range
        last -- the upper index defining the range
        isBP -- True if BP sites, False if LR
        
        Efficiency:
        O(n) where n is the length of the DNA sequence.
    
        Excises elements from a dna sequence by replacing them with 0. Determines the appropriate rec site
        to leave behind.
        '''
        #pdb.set_trace()
        if isBP:
            #P then B on fwd strand
            if sequence[first]%4==2:
                sequence[first]=sequence[first]+8
                for i in range(first+1,last+1):
                    sequence[i]=0
            #B then P on rev strand
            else:
                sequence[last]=sequence[last]+8
                for i in range(first,last):
                    sequence[i]=0
        else:
            #R then L on fwd strand
            if sequence[first]%4==2:
                sequence[first]=sequence[first]-8
                for i in range(first+1,last+1):
                    sequence[i]=0
            #L then R on rev strand
            else:
                sequence[last]=sequence[last]-8
                for i in range(first,last):
                    sequence[i]=0
    
    def invert(self, sequence, first, last, isBP):
        '''
        Inverts elements from a dna sequence between 2 indices.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/23/2018
        
        Arguments:
        sequence -- The dna sequence. A list of int codes
        first -- the lower index defining the range
        last -- the upper index defining the range
        isBP -- True if BP sites, False if LR
        
        Efficiency:
        O(n) where n is the length of the DNA sequence.
    
        Inverts elements from a dna sequence between 2 indices. Determines the appropriate rec sites
        to leave behind.
        '''
        if isBP:
            if sequence[first]%8<3:
                #B then P facing
                if sequence[first]%4==2:
                    sequence[first]=sequence[first]+8
                    sequence[last]=sequence[last]+8
                #B then P not facing
                else:
                    sequence[first]=sequence[first]+12
                    sequence[last]=sequence[last]+4
            else:
                #P then B facing
                if sequence[first]%4==2:
                    sequence[first]=sequence[first]+8
                    sequence[last]=sequence[last]+8
                #P then B not facing
                else:
                    sequence[first]=sequence[first]+4
                    sequence[last]=sequence[last]+12
        else:
            if sequence[first]%8<3:
                #L then R facing
                if sequence[first]%4==2:
                    sequence[first]=sequence[first]-8
                    sequence[last]=sequence[last]-8
                #L then R not facing
                else:
                    sequence[first]=sequence[first]-4
                    sequence[last]=sequence[last]-12
            else:
                #R then L facing
                if sequence[first]%4==2:
                    sequence[first]=sequence[first]-8
                    sequence[last]=sequence[last]-8
                #R then L not facing
                else:
                    sequence[first]=sequence[first]-12
                    sequence[last]=sequence[last]-4
        sequence[first+1:last]=sequence[first+1:last][::-1]
        for i in range(first+1,last):
            if sequence[i]!=0:
                if sequence[i]%4>1:
                    sequence[i]=sequence[i]-2
                else:
                    sequence[i]=sequence[i]+2
                    
    def getProbability(self,window, sequence):
        '''
        Determines the probability of a given recombinase event based on the values from file.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/23/2018
        
        Arguments:
        window -- [start index, end index, is BP?, Recombinase ID, RDF active?, is counter?]
        sequence -- The dna sequence. A list of int codes
        
        Efficiency:
        O(1)
    
        Figures out what type of recombinase event will happen for a given window, and looks up the probability
        that this event happens for that recombinase. Returns this probability.
        '''
        eventType=0
        #if not same direction
        if abs(sequence[window[0]]-sequence[window[1]])!=4:
            eventType+=1
        #if LR
        if not window[2]:
            eventType+=2
        #if the RDF is active
        if window[4]:
            eventType+=4
        return self.recEfficiency[window[3]][eventType]
    
    
    def calculateProbability(self, combo, probabilities, conflicts, conflictGroups,
                             conflictGroupLookup, dependencies, numbCombos, recToWindows):
        '''
        Calculates the probabilities associated with a combination of recombinase windows.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 1/8/2019
        
        Arguments:
        combo -- 
        
        Efficiency:
        O(1)
    
        Figures out what type of recombinase event will happen for a given window, and looks up the probability
        that this event happens for that recombinase. Returns this probability.
        '''
        #pdb.set_trace()
        probs=dict()
        for numbCombo in numbCombos:
            combo2=set(combo)
            conflictGroups2=[]
            probabilities2=probabilities.copy()
            for group in conflictGroups:
                conflictGroups2.append(group.copy())
            for rec in numbCombo:
                for window in recToWindows[rec]:
                    if window in conflictGroupLookup:
                        conflictGroups2[conflictGroupLookup[window]].remove(window)
                    if window in combo2:
                        combo2.remove(window)   
                    probabilities2[window]=0
            prob=1.0
            #changes the window IDs of all windows that are not occuring to their opposites.
            for group in conflictGroups2:
                for i in range(len(group)):
                    if group[i] not in combo2:
                        group[i]=-group[i]
            #Accounts for the probability of the windows that have no conflicts (ie are independent)
            for key in probabilities2:
                if key not in conflictGroupLookup:
                    if key in combo:
                        prob=prob*probabilities2[key]
                    else:
                        prob=prob*(1-probabilities2[key])
            #Calculates the probability of each conflict group. Each group is independent from the other groups and from the independent
            #windows, so once we calculate the probability to the group we can simply mulitply the total probability by it.
            for group in conflictGroups2:
                blocked=False
                #checks that the group is not impeded by another window.
                for dependency in dependencies:
                    if (dependency[0] in group or -dependency[0] in group) and dependency[1] in combo2 and dependency[1] not in group:
                        blocked=True
                        break
                if not blocked:
                    prob=prob*self.resolveConflict(group, probabilities2, conflicts)
            probs[numbCombo]=prob
        return probs
    
    def resolveConflict(self, group, probabilities, conflicts):
        '''
        Returns the probability that a group of recombinase events will occur.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/14/2018
        
        Keyword arguments:
        group -- the group of window IDs
        probabilities -- dictionary of probabilities of each event occuring by itself
        conflicts -- dictionary of conflicting window IDs for each window ID
        
        Efficiency:
        O(n!*n^2) where n is the size of nodes.
        
        Returns the probability that a group of recombinase events will occur. Each event is either assigned to be occuring or not occuring.
        Only possible combinations are considered. Calculates the probability by calculating the probability for each ordering of the events,
        using conditional probabilities for each event to calculate these ordering probabilities. Returns the average of the ordering
        probabilities.
        '''
        #pdb.set_trace()
        #Gets all different orderings of the nodes.
        permutations=itertools.permutations(group)
        prob=0.0
        for perm in permutations:
            prob+=self.subResolve(perm, probabilities, conflicts)
        return prob/math.factorial(len(group))
    
    def subResolve(self, perm, probabilities, conflicts):
        #pdb.set_trace()
        prob2=1.0
        otherProbs=[]
        for i, node in enumerate(perm):
            #Looks up the probability for the node. if the node>0, then that means the window is occuring, if less than 0 it is not
            if node>0:
                prob3=probabilities[node]
            else:
                prob3=1.0-probabilities[-node]
                if i+1<len(perm):
                    for node2 in perm[i+1:]:
                        #case where node depends on node2
                        if node2>0 and node2 in conflicts[abs(node)] and abs(node) not in conflicts[node2]:
                            perm2=list(perm)
                            perm2[i]=-perm2[i]
                            otherProbs.append(self.subResolve(perm2, probabilities, conflicts))
                            break
            #for all previous nodes in the ordering
            for node2 in perm[:i]:
                #if the previous node is occuring, and it direclty conflicts the with the current node
                if node2>0 and node2 in conflicts[abs(node)]:
                    #if the current node is occuring then we have a problem, as the groups should be possible and having both the previous node
                    #and the current node occur is impossible
                    if node>0:
                        print("Error")
                        print(perm, conflicts)
                        raise ValueError("Impossible configuration")
                    #if the previous node is occuring then the current node has a 100% chance of not occuring since they conflict.
                    else:
                        prob3=1
                        break
            prob2=prob2*prob3
        for prob in otherProbs:
            prob2+=prob
        return prob2
    
        #takes a list and generates all possible combinations of the elements of the list
    def getCombinations(self, iterable):
        if len(iterable)>1:
            combinations=self.getCombinations(iterable[:-1])
            combinations2=combinations.copy()
            for c in combinations2:
                c2=c.copy()
                c.append(iterable[-1])
                combinations.append(c2)
            return combinations
        elif len(iterable)==1:
            return [[iterable[0]],[]]
        else:
            return [[]]
    
    #returns all combinations that contain exactly 1 item from each bucket in buckets
    def getCombinations2(self, buckets):
        if buckets:
            bucket=buckets.pop()
            combinations=self.getCombinations2(buckets)
            combinations2=[]
            for b in bucket:
                for c in combinations:
                    c2=c.copy()
                    c2.append(b)
                    combinations2.append(c2)
            return combinations2
        return [[]]
                
    def explore2Way(self, index, conflicts, nodes):
        '''
        Recursive method that finds all nodes in a graph that are connected to a node by 2 way edges.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/14/2018
        
        Keyword arguments:
        index -- the starting index
        conflicts -- dictionary of conflicting window IDs for each window ID
        nodes -- list of interconnected nodes that is being generated by the method
        
        Efficiency:
        O(n) where n is the size of conflicts
        
        Finds all nodes in a graph (conflicts) that are reachable from the current node via 2 way edges. Does not revisit previously
        visited edges.
        '''
        nodes.append(index)
        for node in conflicts[index]:
            try:
                if (node not in nodes) and (index in conflicts[node]):
                    self.explore2Way(node, conflicts, nodes)
            except KeyError:
                pass
    
    
    #def findStates(self, probThreshold, locks, counter, register, currentIndex=0):
    def findStates(self, lock, end2, nextStates, processID, status, finished):
        '''
        Recursively finds the state machine for a genetic circuit.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/14/2018
        
        Keyword arguments:
        probThreshold -- Edges that have a probability below this float will be ignored.
        counter -- list of strings representing a counter dna sequence
        register -- list of strings representing a register dna sequence.
        currentIndex -- The ID of the current circuit configuration.
        
        Efficiency:
        Bad
        
        Recursively finds the state machine for a genetic circuit. Does so by first identifying the active genes,
        then finding pairs of recombinase sites that correspond to those genes, then sorting out conflicts between
        those pairs, then finding the probability that any possible combination of those pairs occurs, then evaluating
        all cases that have probability higher than probThreshold.
        '''
        cycle=False
        #loop is broken when all of the processes
        times={"setup":0.0, "branch":0.0, "find conflicts": 0.0, "numb combos": 0.0, "prob1":0.0, "neighbors":0.0, "prob2": 0.0, "prob3": 0.0, "prob4": 0.0}
        lastWrite=time.time()
        this_process=psutil.Process(os.getpid())
        memory=[]
        readout=open("sm_printout_"+str(processID)+".txt","w")
        while True:
            try:
                memory.append((time.time(), this_process.memory_info().rss/1073741824))
                #if time.time()-240>lastWrite:
                #    readout=open("sm_printout_"+str(processID)+".txt","w")
                #    readout.write(str(times))
                #    readout.close()
                #    lastWrite+=240
                if len(memory)>200:
                    for x in memory:
                        readout.write(str(x[0])+","+str(x[1])+"\n")
                    readout.flush()
                    lastWrite=time.time()
                    memory=[]
                state=nextStates.get(block=True, timeout=10)
                status[processID]=0
                counter=state[0]
                register=state[1]
                currentIndex=state[2]
                #pdb.set_trace()
                expCassettes=set()
                expRDFs=set()
                numb=[]
                
                ##########################################################################
                ##### Identify Active cassettes and RDFs and checks to see if the   ######
                ##### Cycle arrest gene is being expressed.                         ######
                ##########################################################################  
                time1=time.time()
                cycleArrest=self.getExpressedGenes(counter, expCassettes, expRDFs)
                if register:
                    cycleArrest=self.getExpressedGenes(register, expCassettes, expRDFs) or cycleArrest
                if cycleArrest:
                    lock.acquire()
                    end2.send([[currentIndex, [currentIndex,currentIndex], 1],False])
                    end2.recv()
                    lock.release()
                else:
                    ##########################################################################
                    ##### Identifies genes that are fused to Numb.                      ######
                    ##########################################################################  
                    for cassette in expCassettes:
                        if cassette in self.numbIDs:
                            if self.numbIDs[cassette] not in expCassettes:
                                #keeps track of all numb fused genes
                                numb.append(self.numbIDs[cassette])
                                #adds the gene to the list of expressed genes.
                    for key in numb:
                        expCassettes.add(key)
                    ##########################################################################
                    ##### Find all recombinase windows (dna between a pair of sites that #####
                    ##### can perform an action) in the counter and generate all         #####
                    ##### combinations of them that can occur.                           #####
                    ##########################################################################
                    #all the identified windows
                    windows=dict()
                    #used for assigning IDs to windows
                    count=1
                    #the probability of each window in a vacuum
                    probabilities=dict()
                    #lists of all windows for each recombinase
                    recToWindows={x:[] for x in expCassettes}
                    #finds the windows and updates count
                    count=self.getActiveSites(expCassettes, expRDFs, counter, windows, count, True, recToWindows)
                    deleteWindows=[]
                    #figures out the probability for each window
                    for key in windows:
                        windowProb=self.getProbability(windows[key],counter)
                        if windowProb==0:
                            deleteWindows.append(key)
                        else:
                            probabilities[key]=windowProb
                    for key in deleteWindows:
                        del windows[key]
                    time2=time.time()
                    times["setup"]+=-time1+time2
                    conflicts=dict()
                    dependencies=[]
                    combinations={frozenset():0}
                    #finds all combinations of windows where each window can occur simultaneously
                    buckets, windowToBucket=self.branch(windows, conflicts, counter, dependencies, combinations)
                    ##########################################################################
                    ##### Find all recombinase windows (dna between a pair of sites that #####
                    ##### can perform an action) in the register (if there is one) and   #####
                    ##### generate all combinations of these events and the counter      #####
                    ##### that can occur.                                                #####
                    ##########################################################################  
            
                    if register:
                        windows2=dict()
                        count=self.getActiveSites(expCassettes, expRDFs, register, windows2, count, False, recToWindows)
                        for key in windows2:
                            probabilities[key]=self.getProbability(windows2[key],register)
                        buckets2, windowToBucket2=self.branch(windows2, conflicts, register, dependencies, combinations)
                        windows.update(windows2)
                        buckets=buckets+buckets2
                        windowToBucket.update(windowToBucket2)
                    time1=time.time()
                    times["branch"]+=time1-time2
                    memory.append((time.time(), this_process.memory_info().rss/1073741824))
                    ##########################################################################
                    ##### Make lists of mutually exlusive recombinase windows            #####
                    ##########################################################################  
                            
                    conflictGroups=[]
                    conflictGroupLookup=dict()
                    count=0
                    for conflict in conflicts:
                        if conflict not in conflictGroupLookup:
                            deep=[]
                            self.explore2Way(conflict, conflicts, deep)
                            conflictGroups.append(deep)
                            for window in deep:
                                conflictGroupLookup[window]=count
                            count+=1
                    time2=time.time()
                    times["find conflicts"]+=time2-time1
                    numbCombos=[tuple(x) for x in self.getCombinations(numb)]
                    
                    
                    ##########################################################################
                    ##### Determine the probability of each combination of recombinase   #####
                    ##### windows.                                                       #####
                    ##########################################################################  
                    
                    #generate probabilities of all numb outcomes
                    numbProbs=[]
                    numbOutcomes=[]
                    if numb:
                        for n in numb:
                            numbOutcomes.append(((n,0),(n,1)))
                        numbOutcomes=self.getCombinations2(numbOutcomes)
                        numbKeys=[]
                        for outcome in numbOutcomes:
                            numbProbs.append(self.numbProb(outcome))
                            key=[]
                            for x in outcome:
                                #if numb succeeded.
                                if x[1]==0:
                                    key.append(x[0])
                            numbKeys.append(tuple(key))
                    time1=time.time()
                    times["numb combos"]+=time1-time2
                    emptyTuple=tuple()
                    allProbs=[]
                    foundStates=dict()
                    count=0
                    probComparison=[]
                    combinations2=[]
                    #firstCombo=True
                    probSorting=[]
                    memory.append((time.time(), this_process.memory_info().rss/1073741824))
                    for combo in combinations:
                        lastWrite=self.updateMemoryFile(readout, memory, lastWrite, this_process)
                        #calculate probability for combo, for each numb outcome
                        try:
                            prob=self.calculateProbability(combo, probabilities, conflicts,
                                                                  conflictGroups, conflictGroupLookup,
                                                                  dependencies, numbCombos, recToWindows)
                        except ValueError:
                            decoded=self.decode(register)
                            for x in decoded:
                                print(x)
                            print(windows)
                            print("combo",combo)
                            print("conflictGroups", conflictGroups)
                            print("clookup", conflictGroupLookup)
                            print(self.intToText)
                            raise Exception
                        prob2=0.0
                        if numb:
                            for i, key in enumerate(numbKeys):
                                #prob[key] is the probability of the non numb receiving cell, as any x in key is a numb success
                                #so the non numb receiving cell will not receive x. The other cell will always receive all the
                                #recombinases whether or not numb worked, so prob[emptyTuple] (emptyTuple=no missing recs)
                                #gives its probability. A weighted average for each combo is thus calculated, weighted by the
                                #likelihood of this numb outcome. This is used for comparing combos, not for the actual edge
                                #probabilities
                                prob2+=(prob[key]+prob[emptyTuple])*numbProbs[i]
                        else:
                            #if numb is not part of this system, then there is only one probability so no need to average.
                            prob2=prob[emptyTuple]
                        probSorting.append([prob2, SortObject(prob, combo)])
                        #if firstCombo:
                            #firstCombo=False
                            #rbtree=RBTree([prob2, prob, combo])
                        #else:
                            #rbtree.insert([prob2, prob, combo])
                    
                    #rbList=rbtree.toList()
                    #readout.write(str(rbList)+"\n")
                    probSorting.sort()
                    #for i, element in enumerate(rbList):
                    for i, element in enumerate(probSorting):
                        content=element[1].contents()
                        combinations[content[1]]=i
                        allProbs.append(content[0])
                        combinations2.append(content[1])
                        probComparison.append(element[0])
                    #del rbtree
                    #del rbList
                    del probSorting
                    time2=time.time()
                    times["prob1"]+=time2-time1
                    #Find 1st above threshold
                    
                    #check for case where all are good. Need to check because binary search as written will not find the 0th index
                    if probComparison[0]>=self.probThreshold:
                        first_valid_index=0
                    #binary search
                    else:
                        low=0
                        high=len(probComparison)
                        while low+1!=high:
                            mid=(high+low)//2
                            if probComparison[mid]<self.probThreshold:
                                low=mid
                            else:
                                high=mid
                        first_valid_index=high
                    
                    #find an alternative for each combo. We will only consider alterantives that are above the initial threshold. While this may
                    #occasionally miss out on alternaitves that, with added probability, would clear the threshold, the savings in computational time
                    #make it worth it.
                    
                    neighborLists=[]
                    #for each group of recombinase windows
                    foundItem=frozenset()
                    for combo in combinations:
                        lastWrite=self.updateMemoryFile(readout, memory, lastWrite, this_process)
                        #each combo is a neighbor of itself (but we only care about neighbors above initial thresholding)
                        if combinations[combo]>=first_valid_index:
                            neighbors=[combinations[combo]]
                        else:
                            neighbors=[]
                        #cycle through each bucket, finds largest item in bucket contained in the combo
                        #removes this item, then creates neighbors by replacing it with each other item
                        #in the bucket in turn. If no item from bucket is present, does same but w/o
                        #removing an item.
                        for bucket in buckets:
                            combo2=set(combo)
                            length=len(combo)
                            found=False
                            while length>0:
                                for item in bucket:
                                    if len(item)==length:
                                        inCombo2=False
                                        #checks to see if all windows in the group are in combo
                                        for num in item:
                                            if num in combo2:
                                                inCombo2=True
                                            else:
                                                inCombo2=False
                                                break
                                        #if all the windows are in combo
                                        if inCombo2:
                                            #removes each window from combo
                                            for num in item:
                                                combo2.remove(num)
                                            #appends this pared down combo to neighbors
                                            key=frozenset(combo2)
                                            if combinations[key]>=first_valid_index:
                                                neighbors.append(combinations[key])
                                            found=True
                                            foundItem=item
                                            break
                                #add all neighbors that contain a different item from the bucket
                                if found:
                                    for item in bucket:
                                        if item!=foundItem:
                                            combo3=combo2.copy()
                                            for num in item:
                                                combo3.add(num)
                                            key=frozenset(combo3)
                                            if combinations[key]>=first_valid_index:
                                                neighbors.append(combinations[key])
                                    break
                                else:
                                    length+=-1
                            #if there was no item from this bucket in the combo, then we create neighbors with all
                            #items in the bucket
                            if length==0:
                                for item in bucket:
                                    combo3=combo2.copy()
                                    for num in item:
                                        combo3.add(num)
                                    key=frozenset(combo3)
                                    if combinations[key]>=first_valid_index:
                                        neighbors.append(combinations[key])
                        neighborLists.append(neighbors)
                    memory.append((time.time(), this_process.memory_info().rss/1073741824))
                    time1=time.time()
                    times["neighbors"]+=time1-time2
                    #determine actual probability of all pairs of above-threshold-combos.
                    if len(probComparison)<=1000:
                        probCombos=itertools.combinations_with_replacement(range(len(probComparison)), 2)
                        probList=[]
                        probLookup=dict()
                        extraProb=dict()
                        unassigned=0.0
                        #for each pair of combos
                        for p in probCombos:
                            lastWrite=self.updateMemoryFile(readout, memory, lastWrite, this_process)
                            if p[0]<=p[1]:
                                p0=p[0]
                                p1=p[1]
                            else:
                                p0=p[1]
                                p1=p[0]
                            if (p0, p1) in probLookup:
                                prob=probLookup[(p0,p1)]
                            else:
                                state1=allProbs[p0]
                                state2=allProbs[p1]
                                prob=self.calculateFinalProb(state1, state2, p0, p1, numb, numbProbs, numbOutcomes)
                                
                            if p1<first_valid_index:
                                currentBest=None
                                bestValue=prob
                                for neighbor1 in neighborLists[p0]:
                                    for neighbor2 in neighborLists[p1]:
                                        if neighbor1<=neighbor2 and (neighbor1, neighbor2) in probLookup:
                                            tup=(neighbor1, neighbor2)
                                            prob2=probLookup[tup]
                                            
                                        elif neighbor1>neighbor2 and (neighbor2, neighbor1) in probLookup:
                                            tup=(neighbor2, neighbor1)
                                            prob2=probLookup[tup]
                                        else:
                                            state1=allProbs[neighbor1]
                                            state2=allProbs[neighbor2]
                                            prob2=self.calculateFinalProb(state1, state2, neighbor1, neighbor2, numb, numbProbs, numbOutcomes)
                                            if neighbor1<=neighbor2:
                                                tup=(neighbor1, neighbor2)
                                            else:
                                                tup=(neighbor2, neighbor1)
                                            probLookup[tup]=prob2
                                        if prob2>bestValue:
                                            currentBest=tup
                                            bestValue=prob2
                                if currentBest:
                                    if currentBest in extraProb:
                                        extraProb[currentBest]+=prob
                                    else:
                                        extraProb[currentBest]=prob
                                else:
                                    unassigned+=prob
                            elif p0<first_valid_index:
                                currentBest=None
                                bestValue=prob
                                for neighbor1 in neighborLists[p0]:
                                    if neighbor1<=p1 and (neighbor1, p1) in probLookup:
                                        tup=(neighbor1, p1)
                                        prob2=probLookup[tup]
                                        
                                    elif neighbor1>p1 and (p1, neighbor1) in probLookup:
                                        tup=(p1, neighbor1)
                                        prob2=probLookup[tup]
                                    else:
                                        state1=allProbs[neighbor1]
                                        prob2=self.calculateFinalProb(state1, state2, neighbor1, p1, numb, numbProbs, numbOutcomes)
                                        if neighbor1<=p1:
                                            tup=(neighbor1, p1)
                                        else:
                                            tup=(p1, neighbor1)
                                        probLookup[tup]=prob2
                                    if prob2>bestValue:
                                        currentBest=tup
                                        bestValue=prob2
                                if currentBest:
                                    if currentBest in extraProb:
                                        extraProb[currentBest]+=prob
                                    else:
                                        extraProb[currentBest]=prob
                                else:
                                    unassigned+=prob
                            else:
                                tup=(p0,p1)
                                probList.append([prob, tup])
                                if tup not in probLookup:
                                    probLookup[tup]=prob
                    #inaccurate way if accuracy is too computationally costly
                    else:
                        probCombos=itertools.combinations_with_replacement(range(first_valid_index, len(probComparison)), 2)
                        probList=[]
                        probLookup=dict()
                        extraProb=dict()
                        unassigned=1.0
                        #for each pair of combos
                        for p in probCombos:
                            lastWrite=self.updateMemoryFile(readout, memory, lastWrite, this_process)
                            if p[0]<=p[1]:
                                p0=p[0]
                                p1=p[1]
                            else:
                                p0=p[1]
                                p1=p[0]
                            state1=allProbs[p0]
                            state2=allProbs[p1]
                            prob=self.calculateFinalProb(state1, state2, p0, p1, numb, numbProbs, numbOutcomes)
                            unassigned+=-prob
                            tup=(p0,p1)
                            probList.append([prob, tup])
                            probLookup[tup]=prob
                    memory.append((time.time(), this_process.memory_info().rss/1073741824))
                    factor=1-unassigned
                    time2=time.time()
                    times["prob2"]+=time2-time1
                    for p in probList:
                        #reassign prob from neighbors originally thresholded neighbors.
                        if p[1] in extraProb:
                            p[0]=p[0]+extraProb[p[1]]
                        #reassign cut probability proportionally.
                        p[0]=p[0]/factor
                    
                    probList.sort()
                    #reassigns probability to neighbors
                    unassigned=0.0
                    extraProb=dict()
                    for p in probList:
                        lastWrite=self.updateMemoryFile(readout, memory, lastWrite, this_process)
                        #adds accumulated extra probability to the combo (that's been reassigned from other combos)
                        if p[1] in extraProb:
                            p[0]+=extraProb[p[1]]
                        #if probability of combo falls below threshold.
                        if p[0]<self.probThreshold:
                            currentBest=None
                            bestValue=probLookup[p[1]]
                            for neighbor1 in neighborLists[p[1][0]]:
                                for neighbor2 in neighborLists[p[1][1]]:
                                    if neighbor1<=neighbor2:
                                        tup=(neighbor1,neighbor2)
                                    else:
                                        tup=(neighbor2, neighbor1)
                                    prob=probLookup[tup]
                                    if prob>bestValue:
                                        bestValue=prob
                                        currentBest=tup
                            if currentBest:
                                if currentBest in extraProb:
                                    extraProb[currentBest]+=p[0]
                                else:
                                    extraProb[currentBest]=p[0]
                            else:
                                unassigned+=p[0]
                            p.append(False)
                        else:
                            p.append(True)
                    factor=1-unassigned
                    probList.sort()
                    if probList[-1][-1]:
                        found=True
                    time1=time.time()
                    times["prob3"]+=time1-time2
                    #readout.write(str(combinations)+"\n")
                    #readout.write(str(probList)+"\n\n")
                    localEdges=[]
                    localConfigs=dict()
                    #if no combos were above the threshold, cuts combos until the cut probability is enough to make remaining
                    #probs above threshold
                    if not found:
                        index=1
                        found2=False
                        cutProb=0.0
                        while not found2:
                            cutProb+=probList[index-1][0]
                            if (probList[index][0]/(1.0-cutProb))<self.probThreshold:
                                index+=1
                            else:
                                found2=True
                        factor=1-cutProb
                        for p in probList[index:]:
                            p[2]=True
                            p[0]=p[0]/factor
                            for i in p[1]:
                                if i not in foundStates:
                                    comboSet={x: windows[x] for x in combinations2[i]}
                                    counter2=counter.copy()
                                    if register:
                                        register2=register.copy()
                                    else:
                                        register2=None
                                    for key in combinations2[i]:
                                        self.event(comboSet[key], comboSet, counter2, register2)
                                    if register:
                                        key=(tuple(counter2), tuple(register2))
                                    else:
                                        key=(tuple(counter2),)
                                    foundStates[i]=key
                                    localConfigs[key]=0
                    else:
                        for p in probList:
                            if p[2]:
                                #distributes the unaccounted for probability (from thresholded combos with no alternatives)
                                p[0]=p[0]/factor
                                for i in p[1]:
                                    if i not in foundStates:
                                        comboSet={x: windows[x] for x in combinations2[i]}
                                        counter2=counter.copy()
                                        if register:
                                            register2=register.copy()
                                        else:
                                            register2=None
                                        for key in combinations2[i]:
                                            self.event(comboSet[key], comboSet, counter2, register2)
                                        if register:
                                            key=(tuple(counter2), tuple(register2))
                                        else:
                                            key=(tuple(counter2),)
                                        foundStates[i]=key
                                        localConfigs[key]=0
                    time2=time.time()
                    times["prob4"]+=time2-time1
                    #send configs and edges back to main process.
                    lock.acquire()
                    end2.send([localConfigs,True])
                    localConfigs=end2.recv()
                    lock.release()
                    for p in probList:
                        if p[2]:
                            edgeIndices=[]
                            for i in p[1]:
                                edgeIndices.append(localConfigs[foundStates[i]])
                            if len(edgeIndices)==1:
                                edgeIndices.append(edgeIndices[0])
                            localEdges.append([currentIndex, edgeIndices, p[0]])
                    localEdges.append(False)
                    lock.acquire()
                    end2.send(localEdges)
                    end2.recv()
                    lock.release()
                    memory.append((time.time(), this_process.memory_info().rss/1073741824))
                    
                    
            #process exit check
            except queue.Empty:
                status[processID]=1
                if not cycle:
                    allIdle=True
                    for stat in status:
                        if stat==0:
                            allIdle=False
                            break
                    if allIdle:
                        cycle=True
                else:
                    allIdle=True
                    for stat in status:
                        if stat==0:
                            allIdle=False
                            break
                    if allIdle:
                        finished[0]=1
                        readout=open("sm_printout_"+str(processID)+".txt","w")
                        readout.write(str(times))
                        readout.close()
                        try:
                            end2.close()
                            return
                        except:
                            return
    
    def updateMemoryFile(self, readout, memory, lastWrite, process):
        if time.time()-5.0>lastWrite:
            memory.append((time.time(),process.memory_info().rss/1073741824))
            for x in memory:
                readout.write(str(x[0])+","+str(x[1])+"\n")
            readout.flush()
            lastWrite=time.time()
            memory=[]
            return time.time()
        return lastWrite
    
    def calculateFinalProb(self, state1, state2, s1, s2, numb, numbProbs, numbOutcomes):
        emptyTuple=tuple()
        prob=0.0
        if numb:
            for i,outcome in enumerate(numbOutcomes):
                key1=[]
                #for each numb-fused recombinase
                for x in outcome:
                    #if numb succeeded
                    if x[1]==0:
                        key1.append(x[0])
                key1=tuple(key1)
                #key1 is tuple of all recombinases that non-numb-receiving cell did not receive
                #due to success of numb system. emptyTuple means cell got all recombinases.
                if s1==s2:
                    #order does not matter since daugter cells have same final state
                    prob+=state1[key1]*state1[emptyTuple]*numbProbs[i]
                else:
                    #order matters since cell 2 is numb receiving cell, so we calculate both orders and sum
                    prob+=(state1[key1]*state2[emptyTuple]+state2[key1]*state1[emptyTuple])*numbProbs[i]
        else:
            if s1==s2:
                prob=state1[emptyTuple]*state2[emptyTuple]
            else:
                prob=state1[emptyTuple]*state2[emptyTuple]*2
        return prob
    
    #calculates the probability that each combination of numb outcomes will happen. Calculates by assuming that the system has
    #either succeeded or failed, then multiplying that probability by the probability that each outcome has happened given the
    #system wide status.
    def numbProb(self, outcome):
        prob=0.0
        for i in range(2):
            current=self.numbValues[i]
            for x in outcome:
                if x[1]==i:
                    current=current*(self.numbCorrelation[0]+self.numbCorrelation[1]*self.numbValues[i])
                else:
                    current=current*self.numbCorrelation[1]*self.numbValues[x[1]]
            prob+=current
        return prob
                
    def simplify(self, config):
        expCassettes=set()
        genes=[]
        for sequence in config:
            self.getExpressedGenes(sequence, expCassettes, set())
        for cassette in expCassettes:
            name=self.intToText[cassette]
            name=name[:name.index("_")]
            if name.lower() not in self.recombinases:
                numbLoc=name.find("-Numb")
                if (numbLoc!=-1 and (name.lower()[:numbLoc] not in self.recombinases)) or numbLoc==-1:
                    genes.append(name)
        return genes
                
    def searchEdges(self, ids):
        #pdb.set_trace()
        if len(self.allEdges)==0:
            return 0.0
        #pdb.set_trace()
        bottom=0
        top=len(self.allEdges)-1
        while top-bottom>1:
            spot=(top+bottom)//2
            if ids>self.allEdges[spot][:2]:
                bottom=spot
            elif ids<self.allEdges[spot][:2]:
                top=spot
            else:
                return self.allEdges[spot][2]
        if ids==self.allEdges[bottom][:2]:
            return self.allEdges[bottom][2]
        if ids==self.allEdges[top][:2]:
            return self.allEdges[top][2]
        return 0.0
    
    def sortConfigs(self):
        configList=[]
        mapping=dict()
        configs2=dict()
        for config in self.allConfigs:
            configs2[self.allConfigs[config]]=config
        visited=[False]*len(self.allConfigs)
        distances=[]
        for i in range(len(self.allConfigs)):
            distances.append([])
        distances[0].append(0)
        self.exploreConfigs(distances,visited, self.startIndex, 0)
        pairs=[]
        for i in range(len(distances)):
            pairs.append((max(distances[i]),i))
        pairs.sort()
        for i,pair in enumerate(pairs):
            configList.append(configs2[pair[1]])
            mapping[pair[1]]=i
        return configList, mapping
    
    def exploreConfigs(self, distances, visited, index, currentDistance):
        visited[index]=True
        spot=self.searchEdgesByStart(index)
        i=0
        length=len(self.allEdges)
        while spot+i<length and self.allEdges[spot+i][0]==index:
            if self.allEdges[spot+i][1]!=index and not visited[self.allEdges[spot+i][1]]:
                distance=currentDistance+self.allEdges[spot+i][2]
                distances[self.allEdges[spot+i][1]].append(distance)
                self.exploreConfigs(distances, visited.copy(), self.allEdges[spot+i][1], distance)
            i+=1
                
            
    def searchEdgesByStart(self, index):
        #pdb.set_trace()
        if len(self.allEdges)==0:
            return 0
        #pdb.set_trace()
        bottom=0
        top=len(self.allEdges)-1
        while top-bottom>1:
            spot=(top+bottom)//2
            if index>self.allEdges[spot][0]:
                bottom=spot
            elif index<self.allEdges[spot][0]:
                top=spot
            else:
                i=0
                while spot-i-1>-1 and self.allEdges[spot-i-1][0]==index:
                    i+=1
                return spot-i
        if index==self.allEdges[bottom][0]:
            i=0
            while bottom-i-1>-1 and self.allEdges[bottom-i-1][0]==index:
                i+=1
            return bottom-i
        if index==self.allEdges[top][0]:
            i=0
            while top-i-1>-1 and self.allEdges[top-i-1][0]==index:
                i+=1
            return top-i
        raise ValueError("Edge Index not found.")
            
    def run(self, probThreshold, numcores=4):
        self.probThreshold=probThreshold
        lock=Lock()
        end1,end2=Pipe()
        #pdb.set_trace()
        nextStates=Queue()
        l=numcores-1
        finished=Array('i', [0])
        status=Array('i', [0]*l)
        processes=[None]*l
        starttime=time.time()
        this_process=psutil.Process(os.getpid())
        if self.register:
            nextStates.put([self.intCounter.copy(), self.intRegister.copy(), 0])
        else:
            nextStates.put([self.intCounter.copy(), None, 0])
        for i in range(l):
            processes[i]=Process(target=self.findStates,args=(lock, end2, nextStates, i, status, finished))
            processes[i].start()
        while True:
            if time.time()>self.nextTimeStamp:
                self.nextTimeStamp+=30
                print(len(self.allConfigs), len(self.allEdges), nextStates.qsize(), time.time()-starttime, this_process.memory_info().rss)
            isData=end1.poll(5)
            if isData:
                query=end1.recv()
                
                #edges
                if query[-1]==False:
                    end1.send(True)
                    for edge in query[:-1]:
                        self.allEdges.append(edge)
                #configs
                else:
                    response=query[0]
                    for key in response:
                        if key in self.allConfigs:
                            response[key]=self.allConfigs[key]
                        else:
                            response[key]=len(self.allConfigs)
                            self.allConfigs[key]=response[key]
                            if self.register:
                                nextStates.put([list(key[0]), list(key[1]), response[key]])
                            else:
                                nextStates.put([list(key[0]), None, response[key]])
                    end1.send(response)
            elif finished[0]==1:
                end1.close()
                break
        print("found states")
        #open output files
        if "Detailed" in self.outformat:
            self.output=open(self.outFile,"w")
        if "Simple" in self.outformat:
            self.simple=open(self.outFile[:-4]+"_simple.csv", "w")
        self.allEdges.sort()
        if "Detailed" in self.outformat:
            self.output.write("Start Index:\t"+str(self.startIndex)+"\n")
            for config in self.allConfigs:
                self.output.write(str(self.allConfigs[config])+"\t")
                for x in config:
                    self.output.write(str(self.decode(x))+"\t")
                self.output.write("\n")
            self.output.write("\n")
            for edge in self.allEdges:
                self.output.write(str(edge)+"\n")
            self.output.close()
        if "Simple" in self.outformat:
            self.simple.write("Start Index:\t"+str(self.startIndex)+"\n")
            for config in self.allConfigs:
                self.simple.write(str(self.allConfigs[config])+",")
                s=self.simplify(config)
                for x in s:
                    self.simple.write(x+";")
                self.simple.write("\n")
            self.simple.write("\n")
            for edge in self.allEdges:
                self.simple.write(str(edge[0])+",")
                self.simple.write(str(edge[1][0])+",")
                self.simple.write(str(edge[1][1])+",")
                self.simple.write(str(edge[2])+",\n")
            self.simple.close()
        if "Pretty" in self.outformat:
            #configList, mapping=self.sortConfigs()
            #print(mapping)
            configList2=[]
            #for config in configList:
            #    temp=[]
            #    for x in config:
            #        temp.append(self.decode(x))
            #    configList2.append(temp)
            for config in self.allConfigs:
                temp=[]
                for x in config:
                    temp.append(self.decode(x))
                configList2.append(temp)
            edges=[]
            for edge in self.allEdges:
                edges.append([edge[0], edge[1][0], edge[1][1], edge[2]])
            #for edge in self.allEdges:
            #    edge[0]=mapping[edge[0]]
            #    edge[1]=mapping[edge[1]]
            draw.drawStates(configList2, edges, self.pdf)
        print(time.time())

class SortObject():
    def __init__(self, prob, combo):
        self.prob=prob
        self.combo=combo
    
    def __lt__(self, other):
        return True
    
    def __le__(self, other):
        return True
    
    def contents(self):
        return self.prob, self.combo


#red True, black False
#red black tree
class RBTree():
    def __init__(self, firstValue):
        self.top=RBTreeNode(firstValue, False)
    
    def insert(self, value):
        if value[0]>self.top.getValue()[0]:
            nextNode=self.top.getRight()
            if not nextNode:
                newNode=RBTreeNode(value, True, self.top, True)
                self.top.setRight(newNode)
                return
        else:
            nextNode=self.top.getLeft()
            if not nextNode:
                newNode=RBTreeNode(value, True, self.top, False)
                self.top.setLeft(newNode)
                return
        while True:
            if value[0]>nextNode.getValue()[0]:
                nextNode2=nextNode.getRight()
                if nextNode2:
                    nextNode=nextNode2
                else:
                    newNode=RBTreeNode(value, True, nextNode, True)
                    nextNode.setRight(newNode)
                    break
            else:
                nextNode2=nextNode.getLeft()
                if nextNode2:
                    nextNode=nextNode2
                else:
                    newNode=RBTreeNode(value, True, nextNode, False)
                    nextNode.setLeft(newNode)
                    break
        self.repairTree(newNode)
        
    def toList(self):
        l=[]
        self.dfs(self.top, l)
        return l
    
    def dfs(self, node, l):
        left_child=node.getLeft()
        right_child=node.getRight()
        if left_child:
            self.dfs(left_child, l)
        l.append(node.getValue())
        if right_child:
            self.dfs(right_child, l)
        
    def repairTree(self, current):
        #pdb.set_trace()
        parent=current.getParent()
        if parent:
            #red parent
            if parent.getColor():
                grandparent=parent.getParent()
                if parent.getSide():
                    sibling=grandparent.getLeft()
                else:
                    sibling=grandparent.getRight()
                #red uncle (case 3)
                if sibling and sibling.getColor():
                    parent.flipColor()
                    sibling.flipColor()
                    grandparent.setColor(True)
                    self.repairTree(grandparent)
                #black uncle (case 4)
                else:
                    #Downwards: Right, Left
                    if (not current.getSide()) and parent.getSide():
                        self.rotateRight(parent, current, current.getRight())
                        temp=parent
                        parent=current
                        current=temp
                    #Downwards: Left, Right
                    elif (not parent.getSide()) and current.getSide():
                        self.rotateLeft(parent, current, current.getLeft())
                        temp=parent
                        parent=current
                        current=temp
                    if current.getSide():
                        self.rotateLeft(grandparent, parent, parent.getLeft())
                    else:
                        self.rotateRight(grandparent, parent, parent.getRight())
                    parent.setColor(False)
                    grandparent.setColor(True)
            #else case is black parent (case2), no action required
        #is top of tree (case 1)
        else:
            current.setColor(False)
    
    def rotateLeft(self, grandparent, parent, left_child):
        ggp=grandparent.getParent()
        if ggp:
            if grandparent.getSide():
                ggp.setRight(parent)
            else:
                ggp.setLeft(parent)
        else:
            self.top=parent
        parent.setSide(grandparent.getSide())
        parent.setParent(ggp)
        parent.setLeft(grandparent)
        grandparent.setSide(False)
        grandparent.setParent(parent)
        grandparent.setRight(left_child)
        if left_child:
            left_child.setParent(grandparent)
            left_child.setSide(True)
        
    def rotateRight(self, grandparent, parent, right_child):
        ggp=grandparent.getParent()
        if ggp:
            if grandparent.getSide():
                ggp.setRight(parent)
            else:
                ggp.setLeft(parent)
        else:
            self.top=parent
        parent.setSide(grandparent.getSide())
        parent.setParent(ggp)
        parent.setRight(grandparent)
        grandparent.setSide(True)
        grandparent.setParent(parent)
        grandparent.setLeft(right_child)
        if right_child:
            right_child.setParent(grandparent)
            right_child.setSide(False)
        
        
class RBTreeNode():
    def __init__(self, value, color, parent=None, side=False):
        self.value=value
        self.left=None
        self.right=None
        self.color=color
        self.parent=parent
        #False for left, true for right
        self.side=side
    
    def setLeft(self, left=None):
        self.left=left
        
    def setRight(self, right=None):
        self.right=right
    
    def getLeft(self):
        return self.left
    
    def getRight(self):
        return self.right
    
    def getValue(self):
        return self.value
    
    def getColor(self):
        return self.color
    
    def flipColor(self):
        self.color=not self.color
    
    def setColor(self, color):
        self.color=color
        
    def getParent(self):
        return self.parent
    
    def setParent(self, parent=None):
        self.parent=parent
    
    def getSide(self):
        return self.side
    
    def setSide(self, side):
        self.side=side

if __name__=='__main__':
    multiprocessing.set_start_method('fork')
    if 0==1:
        machine=StateMachine(["Simple", "Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachineTest4.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan4count.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Simple", "Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachineTest8.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan8count.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Detailed","Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachineTest4wReg.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan4count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristanUniqueOut4Reg.csv")
        machine.run(.0001)
    if 0==1:
        machine=StateMachine(["Detailed","Simple","Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachineTest2wReg.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2reg.csv")
        machine.run(.0001)
        
    if 0==1:
        machine=StateMachine(["Detailed","Simple","Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachineTest2wRegMulti.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2regmulti.csv")
        machine.run(.0001)
    if 0==1:
        machine=StateMachine(["Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachine4count_2tetexc01.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan4count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2tetbind_exc.csv")
        machine.run(.01)
        
    if 0==1:
        machine=StateMachine(["Simple", "Detailed", "Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachine4count_2tetexc1.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan4count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2tetbind_exc.csv")
        machine.run(.1)
    
    if 0==1:
        machine=StateMachine(["Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachine4count_2tetflip.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan4count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2tetbind_flipreg.csv")
        machine.run(.0001)
    if 0==1:
        machine=StateMachine(["Simple", "Detailed"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachine8count_7tet_1.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan8count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/7tet_example_register.csv")
        machine.run(.1)

    if 0==1:
        machine=StateMachine(["Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachine8count_7tet_002.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan8count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/7tet_example_register.csv")        

        machine.run(.002)
    if 0==1:
        machine=StateMachine(["Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachine8count_7tet_001.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan8count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/7tet_example_register.csv")
        machine.run(.001)
    if 0==1:
        machine=StateMachine(["Detailed","Simple","Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/simple_multi_exc_example.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/simple_multi_exc_example.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter201.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter2.csv")
        machine.run(.01)
        machine=StateMachine(["Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter01.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter_numb01.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter_numb.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Pretty", "Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter_numb201.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter_numb2.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Pretty", "Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter_numb301.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan_example_counter_numb3.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Simple"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachine4count_2tetexc_numb002.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan4count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2tetbind_exc_numb.csv")
        machine.run(.002)
    if 0==1:
        machine=StateMachine(["Pretty"],"examples/meeting_exc_example01.txt",
                             "examples/meeting_exc_example.csv")
        machine.run(.01)
        machine=StateMachine(["Pretty"],"examples/meeting_inv_example01.txt",
                             "examples/meeting_inv_example.csv")
        machine.run(.01)
        machine=StateMachine(["Pretty"],"examples/meeting_rdf_example01.txt",
                             "examples/meeting_rdf_example.csv")
        machine.run(.01)
        machine=StateMachine(["Pretty"],"examples/tristan_example_counter201.txt",
                             "examples/tristan_example_counter2.csv")
        machine.run(.01)
        machine=StateMachine(["Pretty"],"examples/tristan_example_counter01.txt",
                             "examples/tristan_example_counter.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/meeting_multi_rec_example101.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/meeting_multi_rec_example1.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Pretty"],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/meeting_multi_site_example01.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/meeting_multi_site_example.csv")
        machine.run(.01)
    if 0==1:
        machine=StateMachine(["Simple"],"examples/stateMachine8count_7tet_lowerD_numb_early_002_high_mc.txt",
                             "examples/tristan8count.csv",
                             "examples/7tet_lowerD_numb_early.csv")        

        machine.run(.002)
    if 0==1:
        machine=StateMachine(["Simple"],"examples/stateMachine8count_7tet_lowerD_numb_early2_002_max_mc.txt",
                             "examples/tristan8count.csv",
                             "examples/7tet_lowerD_numb_early2.csv","Max")        

        machine.run(.002)
        machine=StateMachine(["Simple"],"examples/stateMachine8count_7tet_lowerD_numb_early2_002_high_mc.txt",
                             "examples/tristan8count.csv",
                             "examples/7tet_lowerD_numb_early2.csv","High")        

        machine.run(.002)
    if 0==1:
        machine=StateMachine(["Simple"],"examples/stateMachine8count_7tet_lowerD_numb_early2_002_mid_mc.txt",
                             "examples/tristan8count.csv",
                             "examples/7tet_lowerD_numb_early2.csv","Mid")        

        machine.run(.002)
        machine=StateMachine(["Simple"],"examples/stateMachine8count_7tet_lowerD_numb_early2_002_low_mc.txt",
                             "examples/tristan8count.csv",
                             "examples/7tet_lowerD_numb_early2.csv","Low")        

        machine.run(.002)
    if 0==1:
        machine=StateMachine(["Simple"],"examples/stateMachine_cone_002_max.txt",
                             "examples/Cone_counter.csv",
                             "examples/Cone_Register.csv", "Max")
        machine.run(.002)
    if 0==1:
        machine=StateMachine(["Simple"],"examples/0homo_3hetero_5cell_circuit_Template.txt",
                             "examples/Cone_counter.csv",
                             "examples/0homo_3hetero_5cell_circuit_Template.csv", "Mid")
        machine.run(.002)
    if 0==1:
        machine=StateMachine(["Simple"],"examples/2018_01_18_lim_green_red_state_machine.txt",
                             "examples/Cone_counter.csv",
                             "examples/Lim_Green_Red.csv", "High")
        machine.run(.002, 8)
    if 1==1:
        machine=StateMachine(["Simple"],"examples/2019_06_04_3tets_register_state_machine.txt",
                             "examples/2019_06_04_counter.csv",
                             "examples/2019_06_04_3tets_register.csv", "High")
        machine.run(.002)