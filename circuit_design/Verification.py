import sys
import itertools
import string
import numpy as np

#sys.path.append('../../')
#from python import helper as hp

sys.path.append('../../python/')
import helper as hp

class Verification():
    
    def __init__(self):
        self.elements = ['promoter_U','promoter_D',
                    'cassette_U', 'cassette_D',
                    'terminator_U', 'terminator_D',
                    'reversal_fac_U', 'reversal_fac_D']

        self.BP_sites = ['BP_site_R', 'BP_site_L']
        self.terminators = ['terminator_U', 'terminator_D']
        self.cassettes = ['cassette_U', 'cassette_D']
        self.promoters = ['promoter_U','promoter_D']
        self.reversal_factors = ['reversal_fac_U', 'reversal_fac_D']
    
    def get_indices_promoters(self, sequence, direction):
        '''return indices of all promoters in the sequence'''
        if direction is 'U':
            return [i for i,x in enumerate(sequence) if 'promoter_U' in x]
        elif direction is 'D':
            return [i for i,x in enumerate(sequence) if 'promoter_D' in x]
        else:
            raise ValueError('Direction must be U or D.')
        
    def get_indices_terminator(self, sequence, direction):
        '''return indices of all terminator in the sequence'''
        if direction is 'U': 
            return [i for i,x in enumerate(sequence) if 'terminator_U' in x]
        elif direction is 'D':
            return [i for i,x in enumerate(sequence) if 'terminator_D' in x]
        else:
            raise ValueError('Direction must be U or D.')
            
    def get_indices_and_ID_cassette(self, sequence, indices_promoters, direction):
        '''return indices of cassette in the sequence after the first promoter'''
        if not indices_promoters:
            return []
        elif direction is 'U': 
            return [(i,x[:2]) for i,x in enumerate(sequence) 
                    if i>indices_promoters[0] and 'cassette_U' in x]
        elif direction is 'D':
            return [(i,x[:2]) for i,x in enumerate(sequence) 
                    if i>indices_promoters[0] and 'cassette_D' in x]
        else:
            raise ValueError('Direction must be U or D.')
            
            
    def get_indices_and_ID_reversal_cofactors(self, sequence, direction):
        '''return indices of all reversal_cofac_U in the sequence'''
        if direction is 'U':
            return [(i,x[:2]) for i,x in enumerate(sequence) if 'reversal_fac_U' in x]
        elif direction is 'D':
            return [(i,x[:2]) for i,x in enumerate(sequence) if 'reversal_fac_D' in x]
        else:
            raise ValueError('Direction must be U or D.')
            
    def get_possible_activated_sites(self, sequence, indices_promoters, indices_terminators, direction, verbose):
        '''find all possible sites that could be activated. ie sites between a promoter and a terminator'''
        # strategy: we look for the first promoter in the list,
        # and find where is the closest terminator on its right
        # return: a list of bool corresponding to the sequence.
        # true if it can be activated, false otherwise.
        
        activated_list = [False]*len(sequence)
        if direction is 'U':
            for idx_promoter in indices_promoters:
                try:
                    idx_terminator = min(filter(lambda x: x > idx_promoter,indices_terminators))
                    # +1 and -1 in next line because we don't care about the promoter,
                    # we're interested in what is in between the promoter and the terminator
                    activated_list[idx_promoter+1:idx_terminator] = [True]*(idx_terminator-idx_promoter-1)
                except:
                    activated_list[idx_promoter+1:] = [True]*len(activated_list[idx_promoter+1:])
                    
        elif direction is 'D':
             for idx_promoter in indices_promoters:
                try:
                    idx_terminator = min(filter(lambda x: x > idx_promoter,indices_terminators))
                    # +1 and -1 in next line because we don't care about the promoter,
                    # we're interested in what is in between the promoter and the terminator
                    activated_list[idx_promoter+1:idx_terminator] = [True]*(idx_terminator-idx_promoter-1)
                except:
                    activated_list[idx_promoter+1:] = [True]*len(activated_list[idx_promoter+1:])
                    
        return activated_list
            

    def get_activate_cassettes(self, indices_and_ID_cassette, activated_sites):
        '''find the cassettes that are activated if any'''
        activated_cassettes = []
        for cassette_idx_ID in indices_and_ID_cassette:
            # we just access the idx
            if activated_sites[cassette_idx_ID[0]]:
                activated_cassettes.append(cassette_idx_ID)
        return activated_cassettes
    
    def get_activate_reversal_cofactors(self, indices_and_ID_reversal_cofactors, activated_sites):
        '''find the reversal_factors that are activated if any'''
        activated_reversal_cofactors = []
        for reversal_cofactors_idx_ID in indices_and_ID_reversal_cofactors:
            # we just access the idx
            if activated_sites[reversal_cofactors_idx_ID[0]]:
                activated_reversal_cofactors.append(reversal_cofactors_idx_ID)
        return activated_reversal_cofactors
            
            
    #TODO there must be a better way to do that than with a nested for loop
    def get_reversal_cofactor_LR_site_ON_idx_ID(self, activated_cassettes, activated_reversal_cofactors, sequence):
        '''find which reversal cofactor is ON: ie we check for the status of the cassette,
        the reversal cofactor and the LR site'''            
        matches = []
        for reversal_cofactor in activated_reversal_cofactors:
            for cassette in activated_cassettes:
                if reversal_cofactor[1] == cassette[1]:
                    # we check if we have LR sites for this id
                    matches.append(reversal_cofactor)
                 
        # look if there is also LR site corresponding to the matches id
        idx_id_ON = []
        
        for match in matches:
            on = [[i,x[:2]] for i,x in enumerate(sequence) if match[1] == x[:2] and 'LR' in x]
            # one rs_site might have been removed by an excision
            if len(on) == 2:
                idx_id_ON.append(on)
                                            
        return idx_id_ON        
            
            
    def get_indices_BP_site_of_cassette(self, sequence, activated_cassette, activated_reversal_cofactors):
        '''return indices of the BP_sites corresponding to the cassette which is activated'''
        indices_bp = [i for i,x in enumerate(sequence) if x[0]==activated_cassette[1][0] and 'BP_site' in x]
                
        if indices_bp:
            # we check if the reversal factor is also active, if yes nothing should happen
            for val_re in activated_reversal_cofactors:
                if activated_cassette[1] in val_re[1]:
                    return []
        return indices_bp            
            
            
    def find_action(self, sequence, BP_sites):
        '''find if we have to perform an excision or an inversion'''    
        # if there is only one bp_sites remaining (eg because of an excision)
        if len(BP_sites) != 2:
            return None
        # if they don't have the same orientation
        elif sequence[BP_sites[0]][-1:]!=sequence[BP_sites[1]][-1:]:
            return 'inversion'
        # if they have the same orientation
        else:
            return 'excision'
        
        
    def get_reversal_action(self, seq, reversal_cofactor_LR_site_ON_idx_ID_i):
        '''find the reversal action'''
        # if they don't have the same orientation   
        if seq[reversal_cofactor_LR_site_ON_idx_ID_i[0][0]][-1]!= \
            seq[reversal_cofactor_LR_site_ON_idx_ID_i[1][0]][-1] :
            return 'inversion'
        # if they have the same orientation
        else:
            return 'excision'     
        
    def get_indices_LR_site_of_reversalco(self, sequence, reversal_cofactor_LR_site_ON_idx_ID_i):
        '''return indices of the BP_sites corresponding to the cassette which is activated'''
        if reversal_cofactor_LR_site_ON_idx_ID_i:
            return [i for i,x in enumerate(sequence) if x[:2] \
                    in reversal_cofactor_LR_site_ON_idx_ID_i[1] and 'LR_site' in x]
        else:
            return []
        
    def excise_sequence(self, sequence, idx_rs_sites_cassettes):
        '''change val to be excise by None'''
        # we don't remove them directly because it
        # would create a conflict with indices.         
        up_seq = sequence[:]
        for idx,_ in enumerate(up_seq):
            if idx >= idx_rs_sites_cassettes[0] and idx < idx_rs_sites_cassettes[1]:
                up_seq[idx] = None
            elif idx == idx_rs_sites_cassettes[1]:
                # turn bp site to lr OR lr to bp if reversal action
                # pass: if an excision is within an excision, sites
                # might already be None, thus the call will raise an error
                if up_seq[idx] == None:
                    pass
                elif 'BP' in up_seq[idx]:
                    up_seq[idx] = up_seq[idx][:3] + '_LR_site_' + up_seq[idx][-1]
                elif 'LR' in up_seq[idx]:
                    up_seq[idx] = up_seq[idx][:3] + '_BP_site_' + up_seq[idx][-1]
                
        return up_seq
    
    
    def inverse_sequence(self, sequence, idx_rs_sites_cassettes):
        '''inverse the part between idx1 and idx2, both included'''
        updated_sequence = sequence.copy()
        idx1 = idx_rs_sites_cassettes[0]
        idx2 = idx_rs_sites_cassettes[1]
        
        # change the orientation letter in each element
        for idx,element in enumerate(sequence):
            if idx1 < idx < idx2:
                if element == None:
                    pass
                elif element[-1] == 'U':
                    updated_sequence[idx] = element[:-1] + 'D'
                elif element[-1] == 'D':
                    updated_sequence[idx] = element[:-1] + 'U'
                elif element[-1] == 'R':
                    updated_sequence[idx] = element[:-1] + 'L'
                elif element[-1] == 'L':
                    updated_sequence[idx] = element[:-1] + 'R' 
                    
        # turn bp sites to lr OR lr to bp (for reversal action)
        if updated_sequence[idx1] == None:
            pass
        elif updated_sequence[idx2] == None:
            pass
        elif 'BP' in updated_sequence[idx1]:
            updated_sequence[idx1] = updated_sequence[idx1][:3] + '_LR_site_' + updated_sequence[idx1][-1]
            updated_sequence[idx2] = updated_sequence[idx2][:3] + '_LR_site_' + updated_sequence[idx2][-1]
        elif 'LR' in updated_sequence[idx1]:
            updated_sequence[idx1] = updated_sequence[idx1][:3] + '_BP_site_' + updated_sequence[idx1][-1]
            updated_sequence[idx2] = updated_sequence[idx2][:3] + '_BP_site_' + updated_sequence[idx2][-1]
            
        # reverse the subsequence being inversed
        updated_sequence = updated_sequence[:idx1+1] + hp.reverse(updated_sequence[idx1+1:idx2]) +\
                                                                updated_sequence[idx2:]
    
        return updated_sequence
    
    
    def are_parallel_sequence_created(self, indices_BP_sites_cassettes):
        '''look if BP_sites are overlapping; if yes, parallel sequences might be created'''
        for a, b in itertools.combinations(indices_BP_sites_cassettes, 2):
            if a[0] > b[0] and a[0] < b[1]:
                return True
            elif a[1] > b[0] and a[1] < b[1]:
                return True
        return False
    
    
    def reverse_idx_to_U(self, rev_idx, sequence):
        '''give the idx position on the seq from the reversed seq'''
        return len(sequence) - rev_idx -1

    def inverse_indices_value_on_U(self, sequence, indices_BP_sites_cassettes_D):
        '''give the indices position on the seq from the reversed seq'''
        indices_on_U = []
        for indices in indices_BP_sites_cassettes_D:
            i_0 = self.reverse_idx_to_U(indices[0],sequence)
            i_1 = self.reverse_idx_to_U(indices[1],sequence)
            indices_on_U.append([i_1, i_0])
        return indices_on_U 
    
    def generate_idx_list_for_fork_actions(self, actions):
        '''return a list of all possible permuatoins of actions'''
        # we permute a list of idx, so we can track
        # which actions correspond to which recombinase
        # think of the example where we have two inversion
        # we will lost the info on which inversion is for which
        # recombinase otherwise
        
        indices = list(range(len(actions)))
        return list(itertools.permutations(indices))
    
    def get_back_parallel_work(self, all_possible_actions, all_actions, all_rs_sites):
        forking_actions = []
        forking_re_sites = []
        for indices in all_possible_actions:
            forking_actions.append([all_actions[i] for i in indices])
            forking_re_sites.append([all_rs_sites[i] for i in indices])
        return forking_actions, forking_re_sites
    
    
    def update_indices_after_action(self, seq, indices_action, indices_BP_sites_cassettes):
        idx_0 = indices_action[0]
        idx_1 = indices_action[1]
        indices_zero = list(range(len(seq)))
        indices_after_inverse = indices_zero[:idx_0+1] + hp.reverse(indices_zero[idx_0+1:idx_1]) + indices_zero[idx_1:]
        new_indices = []
        for indices in indices_BP_sites_cassettes:
            pair_indices = []
            for val in indices:
                pair_indices.append(indices_after_inverse[val])
            new_indices.append(sorted(pair_indices))  
        
        return new_indices
    
    def update_sequence(self, sequence, action, indices_BP_sites_cassettes_i, indices_BP_sites_cassettes):
        '''update the sequence according the the action'''
        if action == 'inversion':
            # update all_indices:
            new_seq = self.inverse_sequence(sequence, indices_BP_sites_cassettes_i)
            indices_BP_sites_cassettes_New = self.update_indices_after_action(sequence, 
                                                                     indices_BP_sites_cassettes_i,
                                                                     indices_BP_sites_cassettes)
            return new_seq, indices_BP_sites_cassettes_New
                
            
        elif action == 'excision':
            # we don't update idx here because the
            # excision is done at the end of the reading;
            # None values are placed as placeholder before that
            new_seq = self.excise_sequence(sequence, indices_BP_sites_cassettes_i)
            return new_seq, indices_BP_sites_cassettes
        
    def remove_none(self, seq):
        '''excise move values from seq'''
        return [x for x in seq if x is not None]
    
    def read_sequence(self, sequence, verbose=0, fork_ON=False, external_count=1, id_fork=None):
        '''read a sequence and find to how much it can count'''
        # verbose = 0: no print. verbose =1: comparison of sequences 
        # verbose = 2: debug mode with more information
        
        history = []
        current_seq = []
        forking_work = []
        # changes_history keep track of the changes
        # of all the bp and rs sites with their ids
        # this is the main ouput of the program as it
        # is used by the register to express genes
        changes_history = []
        updated_seq = sequence[:]
        # start counting at 1 to be consistent with
        # the first report example
        # note that external count is used when a forking happens,
        # and the count is reported from the forking point
        count = external_count
        
        history.append(updated_seq[:])
        while current_seq != updated_seq:
            if verbose in [1,2]:
                print('**********************')
                print('**********************')
                print('current count: ', count, ' to ', count+1)
            current_seq = updated_seq
            # reversed sequence to work with reading from right to left (D)
            sequence_reversed = hp.reverse(current_seq)
    
            ################
            # FORWARD PASS #
            ################        
            
            # find locations of elements
            indices_promoters_U = self.get_indices_promoters(current_seq, 'U')
            indices_terminators_U = self.get_indices_terminator(current_seq, 'U')
            indices_and_ID_cassette_U = self.get_indices_and_ID_cassette(current_seq, indices_promoters_U, 'U')
            indices_and_ID_reversal_cofactors_U = self.get_indices_and_ID_reversal_cofactors(current_seq, 'U')
    
            # find activated site - ie sites that can be subject to changes
            activated_sites_U = self.get_possible_activated_sites(current_seq, 
                                                             indices_promoters_U, 
                                                             indices_terminators_U, 'U',
                                                                 verbose)
            activated_cassettes_U = self.get_activate_cassettes(indices_and_ID_cassette_U, 
                                                           activated_sites_U)
            
                        
            activated_reversal_cofactors_U = self.get_activate_reversal_cofactors(indices_and_ID_reversal_cofactors_U,
                                                                             activated_sites_U)

            if verbose==2:
                print('activated_sites_U: ', activated_sites_U)
                print('activated_cassettes_U: ', activated_cassettes_U)
                print('activated_reversal_cofactors_U: ', activated_reversal_cofactors_U)
            
            #################
            # BACKWARD PASS #
            #################
            
            # find locations of elements
            indices_promoters_D = self.get_indices_promoters(sequence_reversed, 'D')
            indices_terminators_D = self.get_indices_terminator(sequence_reversed, 'D')
            indices_and_ID_cassette_D = self.get_indices_and_ID_cassette(sequence_reversed, indices_promoters_D, 'D')
            indices_and_ID_reversal_cofactors_D = self.get_indices_and_ID_reversal_cofactors(sequence_reversed, 'D')
    
            # find activated site - ie sites that can be subject to changes
            activated_sites_D = self.get_possible_activated_sites(sequence_reversed, 
                                                             indices_promoters_D, 
                                                             indices_terminators_D, 'D',
                                                                 verbose)
            activated_cassettes_D = self.get_activate_cassettes(indices_and_ID_cassette_D, 
                                                             activated_sites_D)
            activated_reversal_cofactors_D = self.get_activate_reversal_cofactors(indices_and_ID_reversal_cofactors_D,
                                                                                  activated_sites_D)
            
            if verbose==2:
                print('activated_sites_D: ', activated_sites_D)
                print('activated_cassettes_D :', activated_cassettes_D)
                print('activated_reversal_cofactors_D: ', activated_reversal_cofactors_D)
            
            ####################
            # STOP CONDITION 1 #
            ####################
        
            # At this point, we can check if anything is activated. 
            # if not, we can break the program
            if not activated_cassettes_U and not activated_reversal_cofactors_D \
                and not activated_cassettes_D and not activated_reversal_cofactors_D:
                if verbose in [1,2]:
                    print('*******************')
                    print('end of the program, STOP condition 1 activated')
                break
                
            ################
            # FORWARD PASS #
            ################ 
            
            reversal_cofactor_LR_site_ON_idx_ID_U = self.get_reversal_cofactor_LR_site_ON_idx_ID(
                                                                        activated_cassettes_U, 
                                                                        activated_reversal_cofactors_U + activated_reversal_cofactors_D,
                                                                        current_seq)   
            
            if verbose==2:
                print('reversal_cofactor_LR_site_ON_idx_ID_U: ', reversal_cofactor_LR_site_ON_idx_ID_U)
            
            # get back lr sites            
            indices_LR_sites_reversalco_U = []
            for val in reversal_cofactor_LR_site_ON_idx_ID_U:
                indices = self.get_indices_LR_site_of_reversalco(current_seq, val)
                if len(indices)==2:
                    indices_LR_sites_reversalco_U.append(indices)
                    
            
            # find locations of BP sites corresponding to expressed recombinase
            indices_BP_sites_cassettes_U = []
            for val in activated_cassettes_U:
                indices = self.get_indices_BP_site_of_cassette(current_seq, val, 
                                                               activated_reversal_cofactors_U + activated_reversal_cofactors_D)
                
                if len(indices)==2:
                    indices_BP_sites_cassettes_U.append(indices)
            
            
            ################
            # BACKWARD PASS #
            ################ 
            
            # find reversal cofactors site for activated reversal co.
            reversal_cofactor_LR_site_ON_idx_ID_D = self.get_reversal_cofactor_LR_site_ON_idx_ID(
                                                                        activated_cassettes_D, 
                                                                        activated_reversal_cofactors_D + activated_reversal_cofactors_U, 
                                                                        sequence_reversed)
            # get back lr sites
            indices_LR_sites_reversalco_D = []
            for val in reversal_cofactor_LR_site_ON_idx_ID_D:
                indices = self.get_indices_LR_site_of_reversalco(sequence_reversed, val)
                if len(indices)==2:
                    indices_LR_sites_reversalco_D.append(indices)
                    
            # find locations of BP sites corresponding to expressed recombinase
            indices_BP_sites_cassettes_D = []
            for val in activated_cassettes_D:
                indices = self.get_indices_BP_site_of_cassette(sequence_reversed, val, 
                                                               activated_reversal_cofactors_D + activated_reversal_cofactors_U)
                if len(indices)==2:
                    indices_BP_sites_cassettes_D.append(indices)
                    
            if verbose==2:
                print('indices_BP_sites_cassettes_U: ', indices_BP_sites_cassettes_U)
                print('indices_BP_sites_cassettes_D: ', indices_BP_sites_cassettes_D)
    
            ################
            #   ACTIONS    #
            ################
                    
            idx_actions_None_U = []
            actions_U = []
            for idx,BP_sites_ in enumerate(indices_BP_sites_cassettes_U):
                action = self.find_action(current_seq, BP_sites_)
                if action is None:
                    idx_actions_None_U.append(idx)
                else:
                    actions_U.append(action)
                    
            idx_actions_None_D = []
            actions_D = []
            for idx,BP_sites_ in enumerate(indices_BP_sites_cassettes_D):
                action = self.find_action(sequence_reversed, BP_sites_)
                if action is None:
                    idx_actions_None_D.append(idx)
                else:
                    actions_D.append(action)
                    
            if verbose==2:
                print('actions_U: ', actions_U)
                print('actions_D', actions_D)
                    
            # remove rs sites in bp list
            indices_BP_sites_cassettes_U = [x for i,x in enumerate(indices_BP_sites_cassettes_U)\
                                            if i not in idx_actions_None_U]
            indices_BP_sites_cassettes_D = [x for i,x in enumerate(indices_BP_sites_cassettes_D)\
                                            if i not in idx_actions_None_D]
                
            reversal_actions_U = []
            if reversal_cofactor_LR_site_ON_idx_ID_U:
                for val in reversal_cofactor_LR_site_ON_idx_ID_U:
                    reversal_actions_U.append(self.get_reversal_action(current_seq, val))
            
            reversal_actions_D = []
            if reversal_cofactor_LR_site_ON_idx_ID_D:
                for val in reversal_cofactor_LR_site_ON_idx_ID_D:
                    reversal_actions_D.append(self.get_reversal_action(sequence_reversed, val))
    
            ################################
            #   PUT EVERYTHING TOGETHER    #
            ################################
            
            D_indices_on_U = self.inverse_indices_value_on_U(current_seq, indices_BP_sites_cassettes_D)
            D_indices_reversal_on_U = self.inverse_indices_value_on_U(current_seq, indices_LR_sites_reversalco_D)
    
            ######################################           
            #TODO: this part is ugly. Rewrite that all stuff
            all_rs_sites_ = [[indices_BP_sites_cassettes_U] + [D_indices_on_U] +\
                            [indices_LR_sites_reversalco_U] + [D_indices_reversal_on_U]]
            all_actions_ = [[actions_U] + [actions_D] +\
                            [reversal_actions_U] + [reversal_actions_D]]
                        
            bio_all_rs_sites = []
            bio_all_actions = []
            for sites,actions in zip(all_rs_sites_,all_actions_):
                if sites and actions:   
                    bio_all_rs_sites = bio_all_rs_sites + sites
                    bio_all_actions = bio_all_actions + actions
            
            all_actions = [x for i,x in enumerate(bio_all_actions) if bio_all_rs_sites[i]]
            all_rs_sites = [x for i,x in enumerate(bio_all_rs_sites) if bio_all_actions[i]]
            
            flat_list_actions = [item for sublist in all_actions for item in sublist]
            flat_rs_sites = [item for sublist in all_rs_sites for item in sublist]
    
            all_rs_sites = flat_rs_sites
            all_actions = flat_list_actions
            ###################################### 
               
            
            if not all_rs_sites or not all_actions:
                if verbose in [1,2]:
                    print('*******************')
                    print('end of the program, STOP condition 2 activated')
                break
                       
            ###########################
            #   FORKING AND ACTIONS   #
            ###########################
            
            # find if a fork is necessary due to cross-talk between actions
            fork = self.are_parallel_sequence_created(all_rs_sites)
            if fork:
                # we don't take the first actions of all_possible_actions
                # as it will be the action processed in this verification
                # run
                all_possible_actions = self.generate_idx_list_for_fork_actions(all_actions)
                #forking_informations = self.get_back_parallel_work(all_possible_actions[1:], 
                #                                              all_actions, 
                #                                              all_rs_sites)
                forking_work.append([current_seq[:], count, len(all_possible_actions)-1])
            
            if fork_ON and count==external_count:
                all_actions_ = []
                all_rs_sites_ = []
                for order in all_possible_actions[id_fork]:
                    all_rs_sites_.append(all_rs_sites[order])
                    all_actions_.append(all_actions[order])
                
                all_rs_sites = all_rs_sites_
                
            flat_rs_sites = list(itertools.chain.from_iterable(all_rs_sites))
            changes_history.append([x for i,x in enumerate(current_seq) if i in flat_rs_sites])
            
            # a fork becomes a break conditions; we don't want to consider
            # sequences that have an inherent probabilistic outcome
            if fork:
                updated_seq, all_rs_sites_New = self.update_sequence(updated_seq, all_actions[0], 
                                                                         all_rs_sites[0], all_rs_sites)
                all_rs_sites = all_rs_sites_New
                if verbose == 2:
                    print('End of analysis because of a fork')
                break
            else:
                for i,action in enumerate(all_actions):
                    updated_seq, all_rs_sites_New = self.update_sequence(updated_seq, action, 
                                                                         all_rs_sites[i], all_rs_sites)
                    all_rs_sites = all_rs_sites_New
                    
            if verbose == 2:
                print('######## INFORMATION ########')
                print('actions_U: ', actions_U)
                print('actions_rev_U: ', reversal_actions_U)
                print('actions_D: ', actions_D)
                print('actions_rev_D: ', reversal_actions_D)
                print('indices_BP_sites_cassettes_U: ', indices_BP_sites_cassettes_U)
                print('D_indices_on_U: ', D_indices_on_U)
                print('indices_LR_sites_reversalco_U: ', indices_LR_sites_reversalco_U)
                print('D_indices_reversal_on_U: ', D_indices_reversal_on_U)
                print('Activated_sites_U:')
                hp.print_changes(current_seq, activated_sites_U, comp_ON=False)
                print('Activated_sites_D:')
                hp.print_changes(sequence_reversed, activated_sites_D, comp_ON=False)
                print('######## END INFORMATION ########')
            
            if verbose==1:
                hp.print_changes(current_seq, updated_seq)
            # remove None values
            updated_seq = self.remove_none(updated_seq)
        
            ####################
            # STOP CONDITION 3 #
            ####################
            
            # some sequence can create infinite loops
            # by coming back to previous states. So, we break
            # the program if we go back to a previous state
            if updated_seq in history:
                break
            
            history.append(updated_seq)
            count+=1
            
        # note: for now, the returned count is the len
        # of the history. In other words, the number of
        # different states can be in.
    
        return len(history), history, forking_work, changes_history
     
    
    
    def generate_sequence(self, max_id, length=10):
        '''generate new sequences'''
        # we'll use alphabet letter as ids
        ids = string.ascii_lowercase[1:]
        elements_for_sequence = []
        
        # first, we create a seed sequence
        # ie we need at least one promoter, one cassette, one terminator and one pair of BP_sites
        # with the same up or down orientation for the promoter, the terminator and the cassette,
        # otherwise the sequence might just not be functional
        id_1 = 'a_'
        unique_id = '0_'
        BP_sites_1 = np.random.choice(self.BP_sites)
        BP_sites_2 = np.random.choice(self.BP_sites)
        elements_for_sequence.append(id_1 + BP_sites_1)
        elements_for_sequence.append(id_1 + BP_sites_2)
        
        orientation = np.random.choice(['_U', '_D'])
        new_cassette = id_1 + 'cassette' + orientation
        new_promoter = unique_id + 'promoter' + orientation
        new_terminator = unique_id + 'terminator' + orientation
            
        elements_for_sequence.append(new_cassette)
        elements_for_sequence.append(new_promoter)
        elements_for_sequence.append(new_terminator)
        
        # we keep track of the current used id by cassette
        # we will need to sample those when we sample a reversal cofactor
        used_ids = ['a']
        ids_for_rs_and_reversal = used_ids[:]
    
        while len(elements_for_sequence) != length:
            # we need at least 3 empty spots to add a combi cassette + BP_sites
            if length-len(elements_for_sequence)>=3 and len(used_ids) < max_id:
                # if all cassette have their reversal and bp sites,
                # we don't sample more of them
                if ids_for_rs_and_reversal:
                    element = np.random.choice(self.elements)
                    # check if the element need a unique id or not
                    if 'promoter' in element or 'terminator' in element:
                        elements_for_sequence.append(unique_id + element)
                    elif 'reversal' in element:
                        id_re = np.random.choice(ids_for_rs_and_reversal)
                        ids_for_rs_and_reversal.remove(id_re)
                        elements_for_sequence.append(id_re + '_' + element)
                       
                    elif 'cassette' in element:
                        id_ = ids[0]
                        elements_for_sequence.append(id_ + '_' + element)
                        elements_for_sequence.append(id_ + '_' + np.random.choice(self.BP_sites))
                        elements_for_sequence.append(id_ + '_' + np.random.choice(self.BP_sites))
                        used_ids.append(id_)
                        ids_for_rs_and_reversal.append(id_)
                        ids = ids.replace(id_, '')
                else:
                    element = np.random.choice(self.promoters + self.terminators + self.cassettes)
                    if 'promoter' in element or 'terminator' in element:
                        elements_for_sequence.append(unique_id + element)  
                    elif 'cassette' in element:
                        id_ = ids[0]
                        elements_for_sequence.append(id_ + '_' + element)
                        elements_for_sequence.append(id_ + '_' + np.random.choice(self.BP_sites))
                        elements_for_sequence.append(id_ + '_' + np.random.choice(self.BP_sites))
                        used_ids.append(id_)
                        ids_for_rs_and_reversal.append(id_)
                        ids = ids.replace(id_, '')
            
            else:
                if ids_for_rs_and_reversal:
                    element = np.random.choice(self.promoters + self.terminators + self.reversal_factors)
                    if 'reversal' in element:
                        id_re = np.random.choice(ids_for_rs_and_reversal)
                        ids_for_rs_and_reversal.remove(id_re)
                        new_element = id_re + '_' + element
                        elements_for_sequence.append(new_element)
                    else:
                        elements_for_sequence.append(unique_id + element)
                else:
                    element = np.random.choice(self.promoters + self.terminators)
                    elements_for_sequence.append(unique_id + element)
                
        np.random.shuffle(elements_for_sequence)
        
        return elements_for_sequence