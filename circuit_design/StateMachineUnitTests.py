import itertools
import csv
import math
import pdb
import struct
import StateMachineDrawings as draw
import time
import unittest
import StateMachine as stateMachine

class TestStateMachine(unittest.TestCase):
    
    def testConstructor(self):
        sm=stateMachine.StateMachine
    
    sm=stateMachine.StateMachine(["Detailed","Simple"],
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/stateMachineTest2wRegMulti.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2count.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/tristan2regmulti.csv")
    
    def testGetExpressedGenes(self):
        expCassettes=set()
        expRDFs=set()
        sequence=
    
    
    
    
    
    def testGetPromoters(self):
        sequence=["a_promoter_F", "b_promoter_R", "abcd"]
        result=self.sm.getPromoters(sequence, "F")
        self.assertEqual(result, [(0, "a")])
        result=self.sm.getPromoters(sequence, "R")
        self.assertEqual(result, [(1, "b")])
        
    def testGetTerminators(self):
        sequence=["a_terminator_F", "b_terminator_R", "abcd"]
        result=self.sm.getTerminators(sequence, "F")
        self.assertEqual(result, [(0, "a")])
        result=self.sm.getTerminators(sequence, "R")
        self.assertEqual(result, [(1, "b")])
        
    def testGetCassettes(self):
        sequence=["a_cassette_F", "b_cassette_R", "abcd"]
        result=self.sm.getCassettes(sequence, "F")
        self.assertEqual(result, [(0, "a")])
        result=self.sm.getCassettes(sequence, "R")
        self.assertEqual(result, [(1, "b")])
        
    def testGetRDFs(self):
        sequence=["a_RDF_F", "b_RDF_R", "abcd"]
        result=self.sm.getRDFs(sequence, "F")
        self.assertEqual(result, [(0, "a")])
        result=self.sm.getRDFs(sequence, "R")
        self.assertEqual(result, [(1, "b")])
    
    def testGetPromotedElements(self):
        sequence=["a_promoter_F", "b_promoter_R", "c_terminator_R", "d_terminator_F","e_promoter_F","f_promoter_F","g_promoter_R", "h_terminator_R"]
        result=self.sm.getPromotedElements(sequence, "F")
        self.assertEqual(result, [False, True, True, False, False, True, True, True])
        result=self.sm.getPromotedElements(sequence[::-1], "R")
        self.assertEqual(result, [True, False,False, True, True, True, False, False][::-1])
        
    def testFilterExpressed(self):
        c=[(0,"a"),(3,"b")]
        p=[False, True, True, True, False]
        result=self.sm.filterExpressed(c,p)
        self.assertEqual(result, [(3,"b")])
    
    def testFindPartnerSites(self):
        siteID="_1"
        sites=[[0,"a","_1","_attB"],[3,"a","_1","_attB"],[4,"a","_1","_attP"],
               [5,"a","_1","_attP"],[6,"a","_0","_attP"],[7,"a","_0","_attB"]]
        index=0
        result=self.sm.findPartnerSites(siteID,sites, index)
        self.assertEqual(result, [[0,"a","_1","_attB"],[3,"a","_1","_attB"],[4,"a","_1","_attP"],
               [5,"a","_1","_attP"]])
        siteID="_0"
        index=4
        result=self.sm.findPartnerSites(siteID,sites, index)
        self.assertEqual(result, [[6,"a","_0","_attP"],[7,"a","_0","_attB"]])
    
    def testGetActiveSites(self):
        expressedCassettes=[(0,"a"), (1,"b")]
        expressedRDFs=[(4,"a")]
        sequence=["a_cassette_R", "b_cassette_R", "c_promoter_R", "c_promoter_F", "a_RDF_F", "d_terminator_F", "a_1_attB_site_F",
                  "a_0_attB_site_F", "a_1_attB_site_R", "a_1_attP_site_R","b_0_attB_site_F", "b_0_attB_site_F", "b_0_attR_site_F",
                  "b_1_attL_site_R", "b_1_attR_site_R", "e_0_attB_site_F", "e_0_attP_site_F"]
        bp,lr=self.sm.getActiveSites(expressedCassettes, expressedRDFs, sequence)
        self.assertTrue([[]])
    
if __name__=='__main__':
    unittest.main()