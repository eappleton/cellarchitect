import math, itertools, pdb
class Test:    
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
    def resolveConflict(self, nodes, probabilities, conflicts):
        '''
        Returns the probability that a group of recombinase events will occur.
        
        Authors: 
        Tristan Daifuku. Church Lab, Harvard Medical School. tristan_daifuku@hms.harvard.edu
        
        Date: 5/14/2018
        
        Keyword arguments:
        nodes -- the group of window IDs
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
        permutations=itertools.permutations(nodes)
        prob=0.0
        for perm in permutations:
            prob+=self.subResolve(perm, probabilities, conflicts)
        return prob/math.factorial(len(nodes))
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
                        raise ValueError("Impossible configuration")
                    #if the previous node is occuring then the current node has a 100% chance of not occuring since they conflict.
                    else:
                        prob3=1
                        break
            prob2=prob2*prob3
        for prob in otherProbs:
            prob2+=prob
        return prob2
#blue 1, purple 2, red 3
probabilities={3:.9, 1:.9, 2:.9}
conflicts={1:set([2,3]), 2:set(), 3:set([1,2])}
conflictGroups=[[1,3]]
conflictGroupLookup={3:0,1:0,}
dependencies=[(1,2)]
numbCombos=[tuple()]
recToWindows={1:3,5:1,9:2}
test=Test()
combo=frozenset([2])
print(test.calculateProbability(combo, probabilities, conflicts, conflictGroups, conflictGroupLookup, dependencies, numbCombos, recToWindows))
combo=frozenset([3])
print(test.calculateProbability(combo, probabilities, conflicts, conflictGroups, conflictGroupLookup, dependencies, numbCombos, recToWindows))
combo=frozenset([1])
print(test.calculateProbability(combo, probabilities, conflicts, conflictGroups, conflictGroupLookup, dependencies, numbCombos, recToWindows))
combo=frozenset()
print(test.calculateProbability(combo, probabilities, conflicts, conflictGroups, conflictGroupLookup, dependencies, numbCombos, recToWindows))