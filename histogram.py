f = open("./50x50wrapx.txt", "r")
eigenvalues = []
for line in f:
    #if float(line)>0:
    eigenvalues.append(float(line))
eigenvalues.sort()
diff=0.0
for i in range(len(eigenvalues)-1):
    if eigenvalues[i+1]-eigenvalues[i]>diff:
        diff=eigenvalues[i+1]-eigenvalues[i]
print diff
f.close()    
