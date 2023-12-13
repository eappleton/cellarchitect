import os, csv
import numpy as np
import matplotlib.pyplot as plt

def plotBindingFreqs(folder):
    files=os.listdir(folder)
    data=dict()
    times=[]
    categories=set()
    for name in files:
        if ".png"!=name[-4:]:
            file=open(folder+name)
            info=list(csv.reader(file))
            timePoint=int(name[6:-4])
            data[timePoint]=dict()
            times.append(timePoint)
            for line in info:
                key=(int(line[0]),int(line[1]))
                categories.add(key)
                data[timePoint][key]=float(line[2])
    times.sort()
    categories=list(categories)
    sums=[]
    lists=[]
    for i in range(len(categories)):
        sums.append([0.0,i])
        lists.append([])
    for time in times:
        for i,key in enumerate(categories):
            if key in data[time]:
                value=data[time][key]
                sums[i][0]+=value
            else:
                value=0.0
            lists[i].append(value)
    sums.sort(reverse=True)
    lists2=[]
    labels=[]
    for sum in sums:
        lists2.append(lists[sum[1]])
        labels.append(str(categories[sum[1]]))
    y=np.vstack(lists2)
    fig, ax=plt.subplots()
    ax.stackplot(times,y,labels=labels)
    ax.legend(loc='upper left', prop={'size': 6})
    plt.savefig(folder+"timeplot.png", dpi=300, format="png")
if __name__=="__main__":
    plotBindingFreqs("simulations_1.0/ape1_6_large_2019-03-07_14-32-48/joint_stats/")