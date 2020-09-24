import numpy as np
import numpy as np
import pylab as pl

lz, lx = 30, 30
N = 20
eigval = [lx, lz]
M = []
data = np.loadtxt('../outputs/string0.out')
for n in range(0,N):
    i=1
    fevalout = open("results/eigenvalues%d.dat" % n,'w')
    fevecout = open("results/eigenvectors%d.dat" % n,'w')

    for layer in range(0,lz*lx):
        Q1 = data[(i-1)*5]
        Q2 = data[((i-1)*5)+1]
        Q3 = data[((i-1)*5)+2]
        Q4 = data[((i-1)*5)+3]
        Q5 = data[((i-1)*5)+4]

        M = np.array([[Q1,Q2,Q3],[Q2,Q4,Q5],[Q3,Q5,-Q1-Q4]])

        eval = np.array(np.linalg.eigvals(M))
        eval, evec = np.linalg.eig(M)
        print(eval[0], evec[:,0])
        i+=1
        fevalout.write(str(eval[0])+"\n")
        fevecout.write(str(evec[0,0])+" "+str(evec[1,0])+"  "+str(evec[2,0])+"\n")

    fevalout.close()
    fevecout.close()
#print(eigval)
