import random, csv, datetime, pdb

def colorConvert(intColor):
    table=['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    outColor='#'
    for c in intColor:
        first=c//16
        second=c%16
        outColor=outColor+table[first]+table[second]
    return outColor

def randomColor(floor=0):
    r=random.randint(0,255)
    g=random.randint(0,255)
    b=random.randint(0,255)
    s=r+g+b
    #while brightness too low
    while s<floor:
        m=min([r,g,b])
        #increase value of lowest channel
        if floor-s+m<255:
            if r==m:
                r=random.randint(floor-s+r,255)
            elif g==m:
                g=random.randint(floor-s+g,255)
            else:
                b=random.randint(floor-s+b,255)
        else:
            if r==m:
                r=255
            elif g==m:
                g=255
            else:
                b=255
        s=r+g+b
    return (r,g,b)

def revComp(seq):
    output=""
    #complements sequence
    for c in seq:
        if c=="A":
            output=output+"T"
        elif c=="T":
            output=output+"A"
        elif c=="C":
            output=output+"G"
        elif c=="G":
            output=output+"C"
        elif c=="a":
            output=output+"t"
        elif c=="t":
            output=output+"a"
        elif c=="c":
            output=output+"g"
        elif c=="g":
            output=output+"c"
        else:
            output=output+'N'
    #reverses sequence
    return output[::-1]

def readParts():
    parts=dict()
    #opens a file which contains part names in column 1, part sequence in column 2, and part annotation in column 3
    reader=csv.reader(open("parts_list.csv"))
    for line in reader:
        #makes fwd and rev versions of each part
        parts[line[0]+"F"]=[line[1],line[2]]
        parts[line[0]+"R"]=[revComp(line[1]), line[2]]
    return parts
        

def makeGB(circuitFile, name):
    #pdb.set_trace()
    colors=dict()
    parts=readParts()
    #read in circuit
    reader=csv.reader(open(circuitFile))
    circuit=[line[0] for line in reader]
    fasta=""
    annotations=[]
    for i, part in enumerate(circuit):
        #add part sequence to whole sequence
        start=len(fasta)+1
        fasta=fasta+parts[part][0]
        end=len(fasta)
        #Colors all related parts the same (eg, all PhiC31 associated parts are colored the same)
        identity=part[:part.find("_")]
        if identity not in colors:
            color=colorConvert(randomColor(255))
            colors[identity]=color
        else:
            color=colors[identity]
        annotations.append([start,end,parts[part][1],part[-1]=="F",color])
        
    length=str(len(fasta))
    months=["JAN",'FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    date=datetime.date.today()
    date=str(date.day)+"-"+months[date.month]+"-"+str(date.year)
    output="LOCUS       "+name[:(27-len(length))]+" "+length+" bp ds-DNA     linear       "+date+"\nDEFINITION  .\nFEATURES             Location/Qualifiers\n"
    for annotation in annotations:
        if annotation[3]:
            output=output+"     misc_feature    "+str(annotation[0])+".."+str(annotation[1])+"\n                     /label=\""+annotation[2]+"\"\n"
        else:
            output=output+"     misc_feature    complement("+str(annotation[0])+".."+str(annotation[1])+")\n                     /label=\""+annotation[2]+"\"\n"
        color=annotation[4]
        output=output+"                     /ApEinfo_revcolor="+color+"\n                     /ApEinfo_fwdcolor="+color+"\n"
    output=output+"ORIGIN\n"
    count=1
    for i in range(0,len(fasta),60):
        strCount=str(count)
        output=output+" "*(9-len(strCount))+strCount+" "+fasta[i:i+60]+"\n"
    output=output+"//"
    outFile=open(circuitFile[:-3]+".gb", "w")
    outFile.write(output)
    outFile.close()
        
            
                
if __name__=='__main__':
    if 1==1:
        makeGB("test_circuit.csv", "Test_Circuit")