# Four types of data file:
# 1. Old unmerged data file, has spaced headers, lacks log columns
# 2. Old merged data file, compact header, lacks log columns
# 3. New unmerged data - nothing to do
# 4. New merged data - nothing to do
# Check seed list is okay for all though
# Adapt old headers to new

import sys

columns=6
headerlength=7
oldcolumns=4

def notempty(x):
    if x!='':
        return True
    else:
        return False

def hashcount(s):
    if len(s) > 0 and s[0]=="#":
        return 1
    else:
        return 0

def fixheader(s):
    s=s.split()
    l=[]
    l.append("# Randomisations:\t"+str(s[1])+";")
    l.append("# X:\t"+str(s[7])+";\tY:\t"+str(s[8])+";\tWrapX:\t"+str(s[9])+";\tWrapY:\t"+str(s[10])+";")
    l.append("# nA:\t"+str(s[3])+";\tnB:\t"+str(s[4])+";\tU:\t"+str(s[2])+";")
    l.append("# Current:"+str(s[5])+";\tGauge:"+str(s[6])+";\tFlux:\t"+str(s[14])+";")
    l.append("# Data points:\t"+str(s[11])+";")
    l.append("# Seeds:\t"+str(s[15])+";")
    l.append("#")
    return l

def splittab(s):
    tabsplit = s.split("\t")
    
    if len(tabsplit) == 1:
        tabsplit = s.split()

    return tabsplit

def mergetabs(l):
    s=""
    for i in range(len(l)):
        s+="\t\t\t"+str(l[i])
    return s


#Open file from command line arguments
filename=sys.argv[1]
f=open(filename,"r")
linelist = [line.strip() for line in f]
f.close()


#Sort out headers
print(sum(map(hashcount, linelist)))
if (sum(map(hashcount, linelist))) == 1:
    #We have old style header, fix
    l1=fixheader(linelist[0])
    linelist=l1+linelist[1:]

#Sort out missing columns (if missing)
l2=linelist[headerlength:]
data=map(splittab, l2)
for i in range(len(data)):
    data[i]=filter(notempty, data[i])
if len(data[0])!=columns:
    #Missing columns
    if len(data[0])!=oldcolumns:
        print "Error: Unknown number of columns: " +str(len(data[0]))
    else:
        #Add dummy columns
        for i in range(len(data)):
            data[i].append("-100")
            data[i].append("-100")

            #Remerge data
            l2=map(mergetabs, data)


linelist=linelist[0:7]+l2

#Rebuild string
data=""
for line in linelist:
    data+=line+"\n"

print data

f=open(filename,"r")
backupdata=f.read()
f.close()
f=open(filename+".bak","w")
f.write(backupdata)
f.close()
f=open(filename,"w")
f.write(data)
f.close()
