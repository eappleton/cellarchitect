from StateMachine import StateMachine
import unittest
class TestStateMachine(unittest.TestCase):
    def testBranch(self):
        #branch(self, windows, conflicts, sequence, dependencies, combinations)
        #Maps each promoter and terminator to negative even numbers, the sites for each recombinase to positive even numbers,
        #the genes that are not RDFs to positive odd numbers, and the RDFs to negative odd numbers. Fwd version of component is
        #always code+2 relative to negative version.
        #[start index, end index, is BP?, Recombinase ID, RDF active?, is counter?]
        tests=[]
        solutions=[]
        sequence=[2,2,6]
        windows={0:[0,2,True,0,False,True], 1:[1,2,True,0,False,True]}
        tests.append([windows,sequence])
        solutions.append({frozenset():0, frozenset({0}):0, frozenset({1}):0})
        sequence=[2,2,6,6]
        windows={0:[0,2,True,0,False,True], 1:[1,2,True,0,False,True], 2:[0,3,True,0,False,True], 3:[1,3,True,0,False,True]}
        tests.append([windows,sequence])
        solutions.append({frozenset():0, frozenset({0}):0, frozenset({1}):0, frozenset({2}):0, frozenset({3}):0})
        sequence=[18,2,2,6,6,22]
        windows={0:[1,3,True,0,False,True], 1:[2,3,True,0,False,True], 2:[1,4,True,0,False,True], 3:[2,4,True,0,False,True], 4:[0,5,True,1,False,True]}
        tests.append([windows,sequence])
        solutions.append({frozenset():0, frozenset({0}):0, frozenset({1}):0, frozenset({2}):0, frozenset({3}):0, frozenset({4}):0})
        sequence=[18,2,2,6,6,22,34,36]
        windows={0:[1,3,True,0,False,True], 1:[2,3,True,0,False,True], 2:[1,4,True,0,False,True], 3:[2,4,True,0,False,True], 4:[0,5,True,1,False,True], 5:[6,7,True,2,False,True]}
        tests.append([windows,sequence])
        solutions.append({frozenset():0, frozenset({0}):0, frozenset({1}):0, frozenset({2}):0, frozenset({3}):0, frozenset({4}):0,
                          frozenset({0,5}):0, frozenset({1,5}):0, frozenset({2,5}):0, frozenset({3,5}):0, frozenset({4,5}):0, frozenset({5}):0})
        for i, solution in enumerate(solutions):
            windows, sequence=tests[i]
            conflicts={}
            dependencies=[]
            combinations={frozenset():0}    
            sm=StateMachine(None,None,None,None,None,True)
            sm.branch(windows, conflicts, sequence, dependencies, combinations)
            print(self.assertEqual(combinations, solution))
unittest.main()