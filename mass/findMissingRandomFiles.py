import os, pdb

#pdb.set_trace()
fileList=os.listdir("simulations_1.0")
circuits=[0]*1365
incomplete=[]
for file in fileList:
    f=str(file)
    index1=f.find("circuit")
    if index1>-1:
        index2=f.find("_",index1)
        number=int(f[index1+7:index2])
        fileList2=os.listdir("simulations_1.0/"+f)
        complete=False
        for file2 in fileList2:
            if "csv.168" in file2:
                complete=True
                circuits[number]+=1
                break
        if not complete:
            incomplete.append(file)
print(circuits)
for i,circuit in enumerate(circuits):
    if circuit<10:
        print(i,circuit)
print(sorted(incomplete))