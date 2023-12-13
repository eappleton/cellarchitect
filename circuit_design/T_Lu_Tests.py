import StateMachine as sm
import itertools

recs=["PhiC31", "Bxb1", "TP901", "int12","PhiRv1", "PhiBT1", "PhiJoe", ]
def test5():
    count=0
    perms=Helper.getPermutations(recs,2)
    for perm in perms:
        print(perm)
        counter=["pEF1alpha_promoter_F"]
        for rec in perm:
            counter.append(rec+"_1_attB_site_F")
            counter.append(rec+"_cassette_F")
            counter.append("sv40_terminator_F")
            counter.append(rec+"_1_attP_site_F")
        counter.append("sv40_terminator_F")
        counterFile=open("/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter5.csv","w")
        for x in counter:
            counterFile.write(x+"\n")
        counterFile.close()
        machine=sm.StateMachine([],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_5state_test.txt",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter5.csv",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_5state.csv")
        machine.run(.0001)
        count+=Helper.count(machine)
    print(count)
def test2():
    count=0
    perms=Helper.getPermutations(recs,1)
    for perm in perms:
        counter=["pEF1alpha_promoter_F"]
        for rec in perm:
            counter.append(rec+"_1_attB_site_F")
            counter.append(rec+"_cassette_F")
            counter.append("sv40_terminator_F")
            counter.append(rec+"_1_attP_site_F")
        counter.append("sv40_terminator_F")
        counterFile=open("/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter2.csv","w")
        for x in counter:
            counterFile.write(x+"\n")
        counterFile.close()
        machine=sm.StateMachine([],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_2state_test.txt",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter2.csv",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_2state.csv")
        machine.run(.0001)
        count+=Helper.count(machine)
    print(count)
def test5():
    count=0
    perms=Helper.getPermutations(recs,2)
    for perm in perms:
        counter=["pEF1alpha_promoter_F"]
        for rec in perm:
            counter.append(rec+"_1_attB_site_F")
            counter.append(rec+"_cassette_F")
            counter.append("sv40_terminator_F")
            counter.append(rec+"_1_attP_site_F")
        counter.append("sv40_terminator_F")
        counterFile=open("/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter5.csv","w")
        for x in counter:
            counterFile.write(x+"\n")
        counterFile.close()
        machine=sm.StateMachine([],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_5state_test.txt",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter5.csv",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_5state.csv")
        machine.run(.0001)
        count+=Helper.count(machine)
    print(count)
def test16():
        count=0
        perms=Helper.getPermutations(recs,3)
        for perm in perms:
            counter=["pEF1alpha_promoter_F"]
            for rec in perm:
                counter.append(rec+"_1_attB_site_F")
                counter.append(rec+"_cassette_F")
                counter.append("sv40_terminator_F")
                counter.append(rec+"_1_attP_site_F")
            counter.append("sv40_terminator_F")
            counterFile=open("/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter16.csv","w")
            for x in counter:
                counterFile.write(x+"\n")
            counterFile.close()
            machine=sm.StateMachine([],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_16state_test.txt",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter16.csv",
                             "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_16state.csv")
            machine.run(.0001)
            count+=Helper.count(machine)
        print(count)
def test65():
    count=0
    perms=Helper.getPermutations(recs,4)
    for perm in perms:
        counter=["pEF1alpha_promoter_F"]
        for rec in perm:
            counter.append(rec+"_1_attB_site_F")
            counter.append(rec+"_cassette_F")
            counter.append("sv40_terminator_F")
            counter.append(rec+"_1_attP_site_F")
        counter.append("sv40_terminator_F")
        counterFile=open("/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter65.csv","w")
        for x in counter:
            counterFile.write(x+"\n")
        counterFile.close()
        machine=sm.StateMachine([],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_65state_test.txt",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter65.csv",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_65state.csv")
        machine.run(.0001)
        count+=Helper.count(machine)
    print(count)
def test326():
    count=0
    perms=Helper.getPermutations(recs,5)
    for perm in perms:
        counter=["pEF1alpha_promoter_F"]
        for rec in perm:
            counter.append(rec+"_1_attB_site_F")
            counter.append(rec+"_cassette_F")
            counter.append("sv40_terminator_F")
            counter.append(rec+"_1_attP_site_F")
        counter.append("sv40_terminator_F")
        counterFile=open("/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter326.csv","w")
        for x in counter:
            counterFile.write(x+"\n")
        counterFile.close()
        machine=sm.StateMachine([],"/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_326state_test.txt",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_Counter326.csv",
                         "/Users/Tristan/Documents/CAD_bio/src/python/circuit_design/examples/T_Lu_326state.csv")
        machine.run(.0001)
        count+=Helper.count(machine)
    print(count)
class Helper():
    def count(machine):
        count=0
        for edge in machine.allEdges:
            if edge[0]==edge[1] and edge[2]==1:
                count+=1
        return count
    def getPermutations(recs,n):
        recs2=recs[:n]
        perms=[]
        for i in range(n+1):
            p=itertools.permutations(recs2,i)
            for perm in p:
                perms.append(perm)
        return perms
                
                    

if __name__=='__main__':
    test326()