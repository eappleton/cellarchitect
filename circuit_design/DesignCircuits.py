import math, pdb, random
terminators=["RBG", "sv40", "BGH"]
terminatorIndex=0
siteIndex=0
def generateStop(direction):
    if random.random()<.35:
        return "Stop2_"+direction
    return "Stop1_"+direction

def designCounter(numDivisions, counterRecs, counterArchitecture):
    global terminators
    global terminatorIndex
    #excision linear counter
    if counterArchitecture=="linear":
        counterSequence=counterRecs
        counter=["CycDep_promoter_F"]
        for i in range(numDivisions-1):
            counter.append(counterSequence[i]+"_1_attB_site_F")
            counter.append(counterSequence[i]+"_cassette_F")
            counter.append(terminators[terminatorIndex]+"_terminator_F")
            terminatorIndex=(terminatorIndex+1)%len(terminators)
            counter.append(counterSequence[i-1]+"_1_attP_site_R")
        counter.append(counterSequence[-1]+"_cassette_F")
        counter.append(terminators[terminatorIndex]+"_terminator_F")
        terminatorIndex=(terminatorIndex+1)%len(terminators)
        return counter
    #binary counter
    else:
        if len(counterRecs)==1:
            counter=["CycDep_promoter_F", counterRecs[0][0][0]+"_cassette_F", terminators[terminatorIndex]+"_terminator_F"]
            terminatorIndex=(terminatorIndex+1)%len(terminators)
            return counter
        counter=["CycDep_promoter_F"]
        for tier in range(len(counterRecs)):
            rec1=counterRecs[tier][0][0]
            counter.insert(1, rec1+"_1_attB_site_F")
            if tier<2:
                counter.append("CycDep_promoter_R")
                counter.append(rec1+"_cassette_F")
                if len(counterRecs[tier])>1:
                    for rec in counterRecs[tier][1:]:
                        counter.append(rec[0]+"_cassette_F")
                counter.append(terminators[terminatorIndex]+"_terminator_F")
                terminatorIndex=(terminatorIndex+1)%len(terminators)
                if tier==1:
                    counter.append(terminators[terminatorIndex]+"_terminator_R")
                    terminatorIndex=(terminatorIndex+1)%len(terminators)
                    counter.append(counterRecs[0][0][0]+"_RDF_R")
            else:
                counter.append(rec1+"_cassette_F")
                if len(counterRecs[tier])>1:
                    for rec in counterRecs[tier][1:]:
                        counter.append(rec[0]+"_cassette_F")
                counter.append(terminators[terminatorIndex]+"_terminator_F")
                terminatorIndex=(terminatorIndex+1)%len(terminators)
                counter.append(counterRecs[tier-1][0][0]+"_RDF_R")
            counter.append(rec1+"_1_attP_site_R")
        counter.append(terminators[terminatorIndex]+"_terminator_F")
        terminatorIndex=(terminatorIndex+1)%len(terminators)
        return counter
    
def designRegister(devTree, numDivisions, counterArchitecture):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/7/2018
    
    Arguments:
    
    devTree - list of lists, with each interior list representing a cell, and providing either daughter cells or expressed proteins.
    numDivisions - The number of divisions required by the developmental tree.
    counterArchitecture - "linear" or "binary"
    
    Purpose: Design an expression register based on a developmental tree. Can be made to be compatible with either a binary or linear counting
    architecture. Uses generic recombinase names. 
    '''
    global siteIndex
    recombinases=[]
    recLookup=dict()   
    #generate generic recombinase names.
    for i in range(len(devTree)):
        recombinases.append("Rec"+str(i))
    expressedStandard=[]
    counterRecs=set()
    counterSequence=[]
    tiers=0
    if counterArchitecture=="linear":
        for i in range(numDivisions):
            #pull the recombinases used in the counter, and add them to the lists of standard recombinases expressed at each count.
            rec=recombinases.pop(0)
            expressedStandard.append([[rec,2]])
            recLookup[rec]=i
            counterRecs.add(expressedStandard[-1][0][0])
            counterSequence.append(expressedStandard[-1][0][0])
    else:
        #find number of recombinases needed for binary counter
        logdiv=math.log2(numDivisions)
        tiers=int(logdiv)
        if logdiv>tiers:
            tiers+=1
        for i in range(tiers):
            #pull the recombinases for the counter. do not add to expressed standard because they are expressed at multiple counts which makes it
            #difficult to use them
            counterSequence.append(recombinases.pop(0))
            counterRecs.add(counterSequence[-1][0])
        for i in range(numDivisions+1):
            expressedStandard.append([])
    
    ##############################################################################################################
    ######                      Build base register                                                         ######
    ##############################################################################################################
    
    register1=descendTree(0, expressedStandard, recombinases, devTree, len(devTree)-1, recLookup)[0]
    checkPromoters(register1,recLookup, expressedStandard, recombinases)
    ##############################################################################################################
    ######                  Build expression register for recombinases needed for base registers            ######
    ##############################################################################################################
    global terminators
    global terminatorIndex
    #checks to see if we need to build the expression register
    hasRecs=False
    for i in range(numDivisions):
        for rec in expressedStandard[i]:
            if rec[0] not in counterRecs:
                hasRecs=True
                break
        if hasRecs:
            break
    recsPerTier=[]
    if hasRecs:
        if counterArchitecture=="linear":
            register2=["CycDep_promoter_F"]
            #we add sites for first and 2nd count regardless of whether there are recombinases at those counts because we want a fwd terminator at count
            #1, which then requires a rev terminator at count 2
            register2.append(counterSequence[0]+"_2_attB_site_F")
            #adds recombinases at 1st count in fwd direction
            for rec in expressedStandard[0]:
                if rec[0] not in counterRecs:
                    register2.append(rec[0]+"_cassette_F")
            #we add terminator regardless of whether there were any recombinases at the 1st count in order to terminate promotion
            register2.append(terminators[terminatorIndex]+"_terminator_F")
            terminatorIndex=(terminatorIndex+1)%len(terminators)
            register2.append(terminators[terminatorIndex]+"_terminator_R")
            terminatorIndex=(terminatorIndex+1)%len(terminators)
            #adds recombinases at 2nd count in rev direction
            for rec in expressedStandard[1]:
                if rec[0] not in counterRecs:
                    register2.append(rec[0]+"_cassette_R")
            #close 1st site
            register2.append(counterSequence[0]+"_2_attP_site_R")
            #add recombinases at remaining counts (can be done in loop because they are not special cases like 1st 2 counts.)
            for i in range(2, numDivisions):
                hasCurrent=False
                for rec in expressedStandard[i]:
                    if rec[0] not in counterRecs:
                        if hasCurrent:
                            register2.append(rec[0]+"_cassette_R")
                        else:
                            hasCurrent=True
                            register2.insert(1, counterSequence[i-1]+"_2_attB_site_F")
                            register2.append(terminators[terminatorIndex]+"_terminator_R")
                            terminatorIndex=(terminatorIndex+1)%len(terminators)
                            register2.append(rec[0]+"_cassette_R")
                if hasCurrent:
                    register2.append(counterSequence[i-1]+"_2_attP_site_R")
        #binary counter architecture. will use expression register format developed by Evan Appleton
        else:
            for rec in counterSequence:
                recsPerTier.append([[rec, 1]])
            countsRight=[]
            countsLeft=[]
            #split up counts that are on the right of left of the promoter, and figure out what genes are at what counts
            for i in range(numDivisions):
                if i%4==0 or i%4==3:
                    countsRight.append(BinaryCount(i))
                    for rec in expressedStandard[i]:
                        countsRight[-1].addGene(rec[0]+"_cassette_")
                else:
                    countsLeft.append(BinaryCount(i))
                    for rec in expressedStandard[i]:
                        countsLeft[-1].addGene(rec[0]+"_cassette_")
            #central promoter that will flip back and forth every other count. we use stand in sites of a given tier, that will be assigned
            #to actual recombinases later. sites are given same site index so they can be matched
            register2=["tier0_"+str(siteIndex)+"_attB_site_F", "CycDep_promoter_F", "tier0_"+str(siteIndex)+"_attP_site_R"]
            #increase siteIndex so that each site pair has unique identifier so they can be matched.
            siteIndex+=1
            #get right side register
            register2=register2+descendBReg(tiers-1, countsRight, numDivisions)
            #get left side register, and flip because the flipSequence method returns the fwd version
            register2=flipSequence(descendBReg(tiers-1, countsLeft, numDivisions))+register2
            #keep track of site names that have been assigned to each siteIndex (so that they can be assigned to other half of site pair)
            siteLookup=dict()
            #put in actual terminator and site names for stand in names
            for i in range(len(register2)):
                item=register2[i]
                if "_terminator_" in item:
                    register2[i]=terminators[terminatorIndex]+item
                    terminatorIndex=(terminatorIndex+1)%len(terminators)
                elif "_site_" in item:
                    loc1=item.find("_")
                    loc2=item.find("_att")
                    site_tier=int(item[4:loc1])
                    #unique site identifier, so pairs of sites can be matched
                    site_index=int(item[loc1+1:loc2])
                    if site_index not in siteLookup:
                        if recsPerTier[site_tier][-1][1]==6:
                            recsPerTier[site_tier].append([recombinases.pop(0),1])
                        else:
                            recsPerTier[site_tier][-1][1]+=1
                        siteLookup[site_index]=recsPerTier[site_tier][-1][0]+"_"+str(recsPerTier[site_tier][-1][1])
                    register2[i]=siteLookup[site_index]+item[loc2:]
    #if no recombinases need to be expressed
    else:
        register2=[]
    if counterArchitecture=="linear":
        return register1+register2, counterSequence
    else:
        return register1+register2, recsPerTier
     
def descendBReg(tier, window, numDivisions):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/14/2018
    
    Arguments:
    
    tier - the current tier of recombinase in the binary counter. Lower tier recombinases are the ones that flip more frequently.
    window - List of Count objects that we are considering.
    numDivisions - The total number of divisions of the cell system.
    
    Purpose - Recursively build the binary recombinase expression register. Format of register developped by Evan Appleton, with counts
    0, 3 (reversed), 4, 7 (reversed), ... to the right of the cycle dependent promoter, and counts 1 (reversed), 2, ... to the left. This
    method only makes one side of the register.
    '''
    global siteIndex
    #pdb.set_trace()
    #if there is only one count in the window. this can only happen if this is the final count on this side of the promoter.
    if len(window)==1:
        forReturn=[]
        #The count will always be the 0th, 4th, 8th,... or 1st, 5th, 9th,... count, as the other counts cannot appear by themselves (they would
        #always be paired with one of the counts listed previously). Therefore, we can eliminate the sites immediately surrounding this count,
        #ie the tier 1 sites (the recombinase that flips every 4 counts). However this means that we need to flip every count that is not a multiple,
        #or 1+ a multiple of 8 (I believe this pattern holds) in order to have the counts end up in the correct orientation when they should be 
        #expressed.
        if window[0].count%8<2:
            for gene in window[0].genes:
                forReturn.append(gene+"F")
            forReturn.append("_terminator_F")
        else:
            forReturn.append("_terminator_R")
            for gene in window[0].genes:
                forReturn.append(gene+"R")
        return forReturn
    else:
        #checks if any of the counts in the window have genes. If there are no genes, we can simply put a fwd and rev terminator in with no sites.
        hasGenes=False
        for i in range(len(window)):
            if window[i].hasGenes():
                hasGenes=True
                break
        if hasGenes:
            forReturn=["tier"+str(tier)+"_"+str(siteIndex)+"_attB_site_F"]
            #if we are at the lowest tier
            if tier==1:
                for gene in window[0].genes:
                    forReturn.append(gene+"F")
                forReturn.append("_terminator_F")
                forReturn.append("_terminator_R")
                for gene in window[1].genes:
                    forReturn.append(gene+"R")
            else:
                midPoint=2**(tier-1)
                #go down a tier.
                if len(window)>midPoint:
                    forReturn=forReturn+descendBReg(tier-1, window[0:midPoint], numDivisions)+descendBReg(tier-1, window[midPoint:], numDivisions)
                #deals with case where window is shorter than standard length for this tier (ie it is the last window), and is so short that it
                #cannot be split in 2.
                else:
                    forReturn=forReturn+descendBReg(tier-1, window)
            forReturn.append("tier"+str(tier)+"_"+str(siteIndex)+"_attP_site_R")
            siteIndex+=1
            return forReturn
        else:
            return ["_terminator_F", "_terminator_R"]

def flipSequence(sequence):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/17/2018
    
    Arguments:
    
    sequence - list of strings representing circuit components
    
    Purpose: Returns the reverse complement of the sequence.
    '''
    opposite={"F":"R", "R":"F"}
    sequence=sequence[::-1]
    for i in range(len(sequence)):
        sequence[i]=sequence[i][:-1]+opposite[sequence[i][-1]]
    return sequence
              
def descendTree(currentDivision,expressedStandard, recombinases, devTree, index, recLookup, promoters=dict()):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/17/2018
    
    Arguments:
    
    currentDivision - The current count.
    expressedStandard - non numb-tagged recombinases expressed at each count.
    recombinases - list of unused recombinase names
    devTree - list representing developmental tree, from the newDevTree module.
    index - location in dev tree
    promoters - dictionary holding promoter regions of register for expression at specific counts
    
    Purpose: recursively explore the developmental tree, building the surface protein expression register.
    '''
    #pdb.set_trace()
    global terminators
    global terminatorIndex
    #No Division
    if devTree[index][0]==-1:
        previousDivision=currentDivision-1
        #if a promoter sequence has not yet been assigned to this division
        if previousDivision not in promoters:
            #each recombinase can have at most 6 site pairs
            if len(expressedStandard[previousDivision])==0 or expressedStandard[previousDivision][-1][1]==6:
                rec=recombinases.pop(0)
                recLookup[rec]=previousDivision
                expressedStandard[previousDivision].append([rec,1])
            else:
                rec=expressedStandard[previousDivision][-1][0]
                expressedStandard[previousDivision][-1][1]+=1
            promoters[previousDivision]=[terminators[terminatorIndex]+"_terminator_R", rec+"_"+str(expressedStandard[previousDivision][-1][1])+"_attB_site_F",
                   "pEF1alpha_promoter_R", terminators[(terminatorIndex+1)%len(terminators)]+"_terminator_F",
                   rec+"_"+str(expressedStandard[previousDivision][-1][1])+"_attP_site_R"]
        genes=promoters[previousDivision].copy()
        for protein in devTree[index][1]:
            genes.append(protein+"_cassette_F")
        genes.append(terminators[terminatorIndex]+"_terminator_F")
        terminatorIndex=(terminatorIndex+1)%len(terminators)
        #we return true, because there is only one promoter currently in the circuit
        return genes, True, 0
    #Asymmetrical Division (using multi recombination, since we don't have a way of doing true asymmetrical division.)
    elif devTree[index][0]==-3:
        #each recombinase can have at most 6 site pairs
        if len(expressedStandard[currentDivision])==0 or expressedStandard[currentDivision][-1][1]==6:
            rec=recombinases.pop(0)
            expressedStandard[currentDivision].append([rec,1])
        else:
            rec=expressedStandard[currentDivision][-1][0]
            expressedStandard[currentDivision][-1][1]+=1
        #go down a level in the tree
        components1, singlePromoter1, height1=descendTree(currentDivision+1, expressedStandard, recombinases, devTree, devTree[index][1], recLookup, promoters)
        components2, singlePromoter2, height2=descendTree(currentDivision+1, expressedStandard, recombinases, devTree, devTree[index][2], recLookup, promoters)
        #if the two daughter cells are both using the same promoter and only have 1 promoter
        if components1[1]==components2[1] and singlePromoter1 and singlePromoter2:
            shared=set()
            count=0
            #we take the promoter region
            forReturn=components1[:5]
            attBloc2=len(components2)
            #we look for the first attB site after the promoter region in the second sequence
            i=5
            for component in components2[5:]:
                if "attB_site" in component:
                    attBloc2=i
                    break
                i+=1
            #we take out the genes that are after the promoter but before the found attB site (if one was found) in the 2nd sequence to match to
            #genes in the first sequence. Anything after the attB site cannot be matched, as it does not apply to all daughter cells of the 2nd cell.
            subsetComponents2=components2[5:attBloc2]
            #match genes in the first sequence to genes in the 2nd
            for component in components1[5:]:
                #we stop when we hit the first attB site in the first sequence for the same reason as described above.
                if "attB_site" in component:
                    #count+2 to make sure we trigger the condition that makes the first sequence be placed 1st in the final circuit
                    count+=2
                    break
                #we add matched genes to the circuit
                if component in subsetComponents2 and "terminator" not in component:
                    forReturn.append(component)
                    shared.add(component)
                #we count the number of elements in the 1st sequence that are not in the 2nd.
                else:
                    count+=1
            #we add sites for multi-recombination.
            forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attB_site_F")
            forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attP_site_F")
            if count>1:
                #we add all components in the 1st sequence that were not in the second.
                for component in components1[5:]:
                    if component not in shared:
                        forReturn.append(component)
                #we add p site for asymmetrical division
                forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attP_site_F")
                count2=0
                #we add all the components in components2, tracking the number of components we add
                for component in components2[5:]:
                    if component not in shared:
                        forReturn.append(component)
                        count2+=1
                #if we only added 1 component, (which would be the terminator), we get rid of the terminator before the attP site.
                if count2==1:
                    forReturn.pop(-3)
            #if there are no components in the first sequence (except for the terminator) that are not in the 2nd sequence, we put the first sequence
            #second becasue this lets us save on a terminator (similarly as 2 lines above)
            else:
                for component in components2[5:-1]:
                    if component not in shared:
                        forReturn.append(component)
                forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attP_site_F")
                for component in components1[5:]:
                    if component not in shared:
                        forReturn.append(component)
            return forReturn, True, height1+1
        #if the two daughter cells are using different promoters (ie, they hit cycleArrest at different counts) or using multiple promoters (which
        #would mean that they have descendants that are hitting cycleArrest at different counts). if this is the case, we cannot match any genes, so
        #we just put the sequences one after another, with asymmetrical division sites around the first one. Note that the shorter part of the tree must
        #come first.
        else:
            if height1<height2:
                forReturn=[rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attB_site_F"]
                forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attP_site_F")
                forReturn=forReturn+components1
                forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attP_site_F")
                forReturn=forReturn+components2
                return forReturn, False, height2+1
            else:
                forReturn=[rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attB_site_F"]
                forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attP_site_F")
                forReturn=forReturn+components2
                forReturn.append(rec+"_"+str(expressedStandard[currentDivision][-1][1])+"_attP_site_F")
                forReturn=forReturn+components1
                return forReturn, False, height1+1
    #Standard Division
    else:
        data=descendTree(currentDivision+1, expressedStandard, recombinases, devTree, devTree[index][1], promoters)
        return data[0], data[1], data[2]+1
    
def checkPromoters(sequence,recLookup, expressedStandard, recombinases, promoter="pEF1alpha"):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 5/2/2019
    
    Arguments:
    
    sequence - the register sequence, stored as a list of strings.
    recLookup - dictionary storing the count associated with each recombinase
    expressedStandard - non numb-tagged recombinases expressed at each count.
    recombinases - list of unused recombinase names
    promoter - the promoter to search for
    
    Purpose: Fixes issues associated with having multiple constitutive promoters active at a given count.
    '''
    foundPromoterSites=[]
    #identifies all instances of the promoter
    prom_string=promoter+"_promoter_R"
    for i,x in enumerate(sequence):
        if prom_string==x:
            foundPromoterSites.append([sequence[i-1],i-1])
    #pdb.set_trace()
    #looks for promoters at the same count
    for i in range(len(foundPromoterSites)-1):
        site=foundPromoterSites[i][0]
        for j in range(i+1,len(foundPromoterSites)):
            site2=foundPromoterSites[j][0]
            #if the promoters are active at the same count
            if site==site2:
                #change sites around the 2nd promoter so that they are different from the sites around the 1st
                index=foundPromoterSites[j][1]
                division=recLookup[site[:site.index("_")]]
                if len(expressedStandard[division])==0 or expressedStandard[division][-1][1]==6:
                    rec=recombinases.pop(0)
                    recLookup[rec]=division
                    expressedStandard[division].append([rec,1])
                else:
                    rec=expressedStandard[division][-1][0]
                    expressedStandard[division][-1][1]+=1
                foundPromoterSites[j][0]=rec+"_"+str(expressedStandard[division][-1][1])+"_attB_site_F"
                sequence[index]=foundPromoterSites[j][0]
                sequence[index+3]=rec+"_"+str(expressedStandard[division][-1][1])+"_attP_site_R"
                #place excision sites around the second promoter that are active if the first promoter is not excised.
                addSite1=index-1
                while "attP" not in sequence[addSite1]:
                    addSite1+=-1
                if len(expressedStandard[division-1])==0 or expressedStandard[division-1][-1][1]==6:
                    rec=recombinases.pop(0)
                    recLookup[rec]=division-1
                    expressedStandard[division-1].append([rec,1])
                else:
                    rec=expressedStandard[division-1][-1][0]
                    expressedStandard[division-1][-1][1]+=1
                sequence.insert(index+4, rec+"_"+str(expressedStandard[division-1][-1][1])+"_attP_site_F")
                sequence.insert(addSite1, rec+"_"+str(expressedStandard[division-1][-1][1])+"_attB_site_F")
                foundPromoterSites[j][1]+=1
                if j+1<len(foundPromoterSites):
                    for k in range(j+1, len(foundPromoterSites)):
                        foundPromoterSites[k][1]+=2
                
                    
class BinaryCount():
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/17/2018
    
    Purpose: Simple object that holds a list of genes and a count
    '''
    def __init__(self, count):
        self.genes=[]
        self.count=count
    
    def addGene(self, gene):
        self.genes.append(gene)
        
    def hasGenes(self):
        if self.genes:
            return True
        return False
    
    def __str__(self):
        return str(self.count)+": "+str(self.genes)

def componentEqual(component1, component2):
    '''
    Author: Tristan Daifuku, Harvard Medical School
    Date: 12/17/2018
    
    Arguments:
    
    component1 - string representation of the 1st component.
    component2 - string representation of the 2nd component.
    
    Purpose: Compare 2 circuit components. If they are the same, returns true, or if they are both terminators on the same strand, returns true.
    '''
    if component1==component2:
        return True
    if "terminator" in component1 and "terminator" in component2:
        if component1[-1]==component2[-1]:
            return True
    return False