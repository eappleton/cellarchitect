import os
import shutil

names=set()
fileList=os.listdir("examples/Random/Batch1")
for file in fileList:
    f=str(file)
    index1=f.find("circuit")
    index2=f.find("_",index1)
    names.add(f[index1:index2])
fileList=os.listdir("examples/Random/Batch2")
for file in fileList:
    f=str(file)
    index1=f.find("circuit")
    index2=f.find("_",index1)
    names.add(f[index1:index2])
fileList=os.listdir("simulations_1.0")
count=0
for file in fileList:
    f=str(file)
    #2print(f)
    index1=f.find("circuit")
    index2=f.find("_",index1)
    if f[index1:index2] in names:
        count+=1
        if count%100==0:
            print(count)
        filelist2=os.listdir("simulations_1.0/"+f)
        incomplete=True
        for f2 in filelist2:
            if "csv.168" in str(f2):
                incomplete=False
                break
        if incomplete:
            shutil.copytree("simulations_1.0/"+f,"simulations_1.0/Transfer1/Incomplete/"+f)
        else:
            shutil.copytree("simulations_1.0/"+f, "simulations_1.0/Transfer1/"+f)