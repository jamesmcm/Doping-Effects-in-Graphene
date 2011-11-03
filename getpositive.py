f = open("./50x50wrapx.txt", "r")
fw= open("./50x50wrapxpositivesort.txt", "w")
eigenvalues = []
for line in f:
    if float(line)>0:
        eigenvalues.append(float(line))
eigenvalues.sort()
for i in range(len(eigenvalues)):
    fw.write(str(i) + "\t" + str(eigenvalues[i]) +"\n")
fw.close()
f.close()
    
    
