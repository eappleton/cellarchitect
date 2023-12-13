import matplotlib.pyplot as plt
plt.figure(dpi=300, figsize=(8,5))
for i in range(7):
    f=open("sm_printout_"+str(i)+".txt")
    x=[]
    values=[]
    for line in f:
        a=line.split(",")
        x.append(float(a[0]))
        values.append(float(a[1]))
    
    plt.scatter(x, values, marker=".", s=2)