import time, random
group1=[]
group2=[]
set1=set()
set2=set()
for i in range(100):
    set1.add(random.randint(0,100))
    set2.add(random.randint(100,200))
for i in range(100):
    group1.append(set1)
    group2.append(set2)
runs=10000
start=time.time()
for i in range(runs):
    forRemoval=[]
    count=0
    for x in group1[0]:
        forRemoval.append(x)
        count+=1
        if count==20:
            break
    for x in forRemoval:
        group1[0].remove(x)
    forRemoval=[]
    count=0
    for x in group2[0]:
        forRemoval.append(x)
        count+=1
        if count==20:
            break
    for x in forRemoval:
        group2[0].remove(x)
    for i in range(20):
        group1[0].add(random.randint(0,100))
        group2[0].add(random.randint(100,200))
    
    c=group1[0].union(group2[0])
    for i in range(100):
        group1[i]=c
        group2[i]=c
end=time.time()
print("union reassign",end-start)
group1=[]
group2=[]
set1=set()
set2=set()
for i in range(100):
    set1.add(random.randint(0,100))
    set2.add(random.randint(100,200))
for i in range(100):
    group1.append(set1)
    group2.append(set2)
start=time.time()
for i in range(runs):
    forRemoval=[]
    count=0
    for x in group1[0]:
        forRemoval.append(x)
        count+=1
        if count==20:
            break
    for x in forRemoval:
        group1[0].remove(x)
    forRemoval=[]
    count=0
    for x in group2[0]:
        forRemoval.append(x)
        count+=1
        if count==20:
            break
    for x in forRemoval:
        group2[0].remove(x)
    for i in range(20):
        group1[0].add(random.randint(0,100))
        group2[0].add(random.randint(100,200))
    
    for x in group1[0]:
        group2[0].add(x)
    for x in group2[0]:
        group1[0].add(x)
end=time.time()
print("double loop",end-start)
