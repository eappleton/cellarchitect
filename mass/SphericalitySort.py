import EigenVectorSimilarity as ev
import os
from statistics import mean

def sortBySphericality(folderPath):
    fileList=os.listdir(folderPath)
    scores=dict()
    scores2=[]
    maxScore=0
    count=0
    for file in fileList:
        count+=1
        if count%100==0:
            print(count)
        network1,points_list_1,points_norm_1,scaled_diameter_1=ev.read_points(folderPath+"/"+file+"/cell_positions.csv.168")
        if not points_list_1:
            scores2.append([-1,str(file)])
        else:
            score=ev.Sphericity(points_list_1)
            if score>maxScore:
                maxScore=score
            scores2.append([score,str(file)])
    for i in range(len(scores2)):
        if scores2[i][0]==-1:
            scores2[i][0]=maxScore+1
        index1=scores2[i][1].find("circuit")
        index2=scores2[i][1].find("_",index1)
        name=scores2[i][1][index1:index2]
        if name not in scores:
            scores[name]=[]
        scores[name].append(scores2[i][0])
    scores2.sort()
    scores3=[]
    for name in scores:
        scores3.append([mean(scores[name]),name])
    scores3.sort()
    file2=open(folderPath+"/IndividualSort.csv","w")
    for score in scores2:
        file2.write(score[1]+","+str(score[0])+",\n")
    file2.close()
    file3=open(folderPath+"/CircuitSort.csv","w")
    for score in scores3:
        file3.write(score[1]+","+str(score[0])+",\n")
    file3.close()
    
if __name__=="__main__":
    sortBySphericality("2018_11_07_Random_Results/Transfer1")