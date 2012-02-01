#a_n+1 =-a_n + lambda c_n
#b_n+1 =-b_n + lambda d_n
#c_n+1 =-E a_n + (E*lambda - 1) c_n
#d_n+1 =-E b_n + (E*lambda - 1) d_n

#T = 2/((a+d)+(c-b)i)
#lambda = E +- 1

#initial values, 1 pair block:
#generate +ve lambda set first:
#From E=-5 to +5 in steps of 0.01

E=-5.0
posltlist=[]
negltlist=[]
n=0
i=1
p=10
while (E<5.0):
    i=1
    lamb=E+1.0
    a=-1.0
    b=lamb
    c=-1*E
    d=(E*lamb)-1.0
    while (i<p):
        atemp=(-1.0*a) + (lamb*c)
        btemp=(-1.0*b) + (lamb*d)
        ctemp=(-1.0*E*a) + ((E*lamb)-1)*c
        dtemp=(-1.0*E*b) + ((E*lamb)-1)*d
        a=atemp
        b=btemp
        c=ctemp
        d=dtemp
        i=i+1
    val=(4.0*(pow(a+d, 2) + pow(c-b, 2)))/(pow((pow(a+d,2)+pow(c-b,2)),2))
    posltlist.append(val)

    i=1
    lamb=E-1.0
    a=-1.0
    b=lamb
    c=-1*E
    d=(E*lamb)-1.0
    while (i<p):
        atemp=(-1.0*a) + (lamb*c)
        btemp=(-1.0*b) + (lamb*d)
        ctemp=(-1.0*E*a) + ((E*lamb)-1)*c
        dtemp=(-1.0*E*b) + ((E*lamb)-1)*d
        a=atemp
        b=btemp
        c=ctemp
        d=dtemp
        i=i+1
    val=(4.0*(pow(a+d, 2) + pow(c-b, 2)))/(pow((pow(a+d,2)+pow(c-b,2)),2))
    negltlist.append(val)

    print(str(E) + "\t" + str(posltlist[n]) + "\t" + str(negltlist[n]))
    E=E+0.01
    n=n+1
    
    
        
