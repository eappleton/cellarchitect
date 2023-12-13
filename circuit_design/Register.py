import sys
import numpy as np
import itertools
import string

#sys.path.append('../../')
#from python import helper as hp

sys.path.append('../../python/')
import helper as hp

sys.path.append('../../python/circuit_design/')
from Verification import Verification


class Register(Verification):
    
    def __init__(self):
        
        self.elements = ['promoter_U','promoter_D',
                         'cassette_U', 'cassette_D',
                         'terminator_U', 'terminator_D']
        
        
    def get_rs_sites(self, sequence, id_):
        '''return indices of rs sites in the sequence'''
        return [x for x in sequence if id_ in x[:2] and 'site' in x]
    
    def get_indices_rs_sites(self, sequence, id_):
        '''return indices of rs sites in the sequence'''
        return [i for i,x in enumerate(sequence) if id_ in x[:2] and 'site' in x]
    
    def get_indices_all_rs_sites(self, sequence, id_, id_int):
        '''return indices of rs sites in the sequence'''
        return [i for i,x in enumerate(sequence) if id_ in x[:2] and 'site' in x and id_int in x]
    
    def get_indices_and_ID_gene(self, sequence, indices_promoters, direction):
        '''return indices of gene in the sequence after the first promoter'''
        if not indices_promoters:
            return []
        elif direction is 'U': 
            return [(i,x[:2]) for i,x in enumerate(sequence) 
                    if i>indices_promoters[0] and 'gene' in x and '_U' in x]
        elif direction is 'D':
            return [(i,x[:2]) for i,x in enumerate(sequence) 
                    if i>indices_promoters[0] and 'gene' in x and '_D' in x]
        else:
            raise ValueError('Direction must be U or D.')
            
            
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

    
    def read_sequence(self, sequence, actions_list, unique_genes, verbose=0):
        '''read a sequence and find to how much it can count'''
        # verbose = 0: no print. 
        # verbose = 1: comparison of sequences
        # verbose = 2: debug mode with more information
        
        history = []
        current_seq = []
        # changes_history keep track of the changes
        # of all the bp and rs sites with their ids
        # this is the main ouput of the program as it
        # is used by the register to express genes
        changes_history = []
        current_seq = sequence[:]
                         
        # init expression matrix
        n_genes = len(unique_genes)
        expression_matrix = np.zeros([len(actions_list)+1, n_genes]).astype(int)
       
    
        # We start by looking for the initial expression of the sequence
        self.initial_expression(current_seq, expression_matrix, unique_genes)
        history.append(current_seq[:])
        
        for iter_,actions in enumerate(actions_list):
            if verbose in [1,2]:
                print('**********************')
                print('Coming change to: ', actions, ' iter: ', iter_)
            
            ################
            # Apply action #
            ################
            
            all_rs_sites = []
            all_actions = []
            
            ids = []
            for action in actions:
                ids.append(action[:2])
            ids = list(set(ids))
            
            rs_sites = []
            for id_ in ids:
                rs_sites.append(self.get_rs_sites(current_seq, id_))
        
            if verbose == 2:
                print('id_: ', id_)
                print('rs_sites: ', rs_sites)
            
            # if the corresponding rs_sites are not in the register,
            # we keep the action part
            if rs_sites[0]:
                all_rs_sites = []
                for i,id_ in enumerate(ids):
                    #for now, we deal with maximum 3 recombinase sites for one recombinase
                    
                    size = len(rs_sites[i])
                    range_ = int(size/2)
                    for j in range(range_):
                        if len(rs_sites[i]) == 2:
                            all_rs_sites.append(self.get_indices_rs_sites(current_seq, id_))
                        else:
                            all_rs_sites.append(self.get_indices_all_rs_sites(current_seq, id_, str(j+1))) 
                    if size>=4:
                        fork = self.are_parallel_sequence_created(all_rs_sites)
                    if size>6:
                        print('Warning: there is more than 3 BP sites pairs for a given recombinase')

                
                # find actions
                for rs_sites in all_rs_sites:
                    action = self.find_action(current_seq, rs_sites)
                    all_actions.append(action)
                
                if verbose == 2:
                    print('all_rs_sites: ', all_rs_sites)
                    print('all_actions: ', all_actions)
                
                updated_seq = current_seq[:]
                for i,action in enumerate(all_actions):
                    updated_seq, all_rs_sites_New = self.update_sequence(updated_seq, action, 
                                                                         all_rs_sites[i], all_rs_sites)
                    all_rs_sites = all_rs_sites_New
                        
                # remove None values
                updated_seq = self.remove_none(updated_seq)
            
            else:
                updated_seq = current_seq[:]
                
            
            #######################
            # Read for expression #
            #######################
               
            ################
            # FORWARD PASS #
            ################  
            
            # find locations of elements
            indices_promoters_U = self.get_indices_promoters(updated_seq, 'U')
            indices_terminators_U = self.get_indices_terminator(updated_seq, 'U')
            indices_and_ID_gene_U = self.get_indices_and_ID_gene(updated_seq,
                                                                 indices_promoters_U, 'U')
    
            # find activated site - ie sites that can be subject to changes
            activated_sites_U = self.get_possible_activated_sites(updated_seq, 
                                                             indices_promoters_U, 
                                                             indices_terminators_U, 'U', False)
            activated_gene_U = self.get_activate_cassettes(indices_and_ID_gene_U, 
                                                           activated_sites_U)
           
            
            #################
            # BACKWARD PASS #
            #################
            
            # reversed sequence to work with reading from right to left (D)
            sequence_reversed = hp.reverse(updated_seq)
            
            # find locations of elements
            indices_promoters_D = self.get_indices_promoters(sequence_reversed, 'D')
            indices_terminators_D = self.get_indices_terminator(sequence_reversed, 'D')
            indices_and_ID_gene_D = self.get_indices_and_ID_gene(sequence_reversed, 
                                                                 indices_promoters_D, 'D')
    
            # find activated site - ie sites that can be subject to changes
            activated_sites_D = self.get_possible_activated_sites(sequence_reversed, 
                                                             indices_promoters_D, 
                                                             indices_terminators_D, 'D', False)
            activated_gene_D = self.get_activate_cassettes(indices_and_ID_gene_D, 
                                                           activated_sites_D)
            if verbose==2:
                print('activated_gene_D: ', activated_gene_D)
                print('activated_gene_U: ', activated_gene_U)
                
            # We add the activated id to the matrix
            all_activated_genes = [activated_gene_U,activated_gene_D]    
                
            for genes in all_activated_genes:
                if genes:
                    all_ids = [x[1][0] for x in genes]
                    for id_ in all_ids:
                        expression_matrix[iter_+1][unique_genes.index(id_)] = 1
            
            
            history.append(updated_seq)
            
            if verbose == 2:
                hp.print_changes(current_seq, updated_seq, comp_ON=True)
            current_seq = updated_seq
    
        return history, expression_matrix
    
    
    
    def initial_expression(self, sequence, expression_matrix, unique_genes):
        
        # find locations of elements
        indices_promoters_U = self.get_indices_promoters(sequence, 'U')
        indices_terminators_U = self.get_indices_terminator(sequence, 'U')
        indices_and_ID_gene_U = self.get_indices_and_ID_gene(sequence,
                                                             indices_promoters_U, 'U')
    
        # find activated site - ie sites that can be subject to changes
        activated_sites_U = self.get_possible_activated_sites(sequence, 
                                                         indices_promoters_U, 
                                                         indices_terminators_U, 'U', False)
        activated_gene_U = self.get_activate_cassettes(indices_and_ID_gene_U, 
                                                       activated_sites_U)
        
        
        #################
        # BACKWARD PASS #
        #################
        
        # reversed sequence to work with reading from right to left (D)
        sequence_reversed = hp.reverse(sequence)
        
        # find locations of elements
        indices_promoters_D = self.get_indices_promoters(sequence_reversed, 'D')
        indices_terminators_D = self.get_indices_terminator(sequence_reversed, 'D')
        indices_and_ID_gene_D = self.get_indices_and_ID_gene(sequence_reversed, 
                                                                 indices_promoters_D, 'D')
    
        # find activated site - ie sites that can be subject to changes
        activated_sites_D = self.get_possible_activated_sites(sequence_reversed, 
                                                         indices_promoters_D, 
                                                         indices_terminators_D, 'D', False)
        activated_gene_D = self.get_activate_cassettes(indices_and_ID_gene_D, 
                                                         activated_sites_D)

        
        # We add the activated id to the matrix
        all_activated_genes = [activated_gene_U,activated_gene_D]
        
        for genes in all_activated_genes:
            if genes:
                id_ = genes[0][1][0]
                expression_matrix[0][unique_genes.index(id_)] = 1
                        
           
