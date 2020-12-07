import numpy as np
import numpy as np
import pylab as pl

lz, lx = 100, 100
ly = 1
N = 50
eigval = [lx, lz]
M = []
data = np.loadtxt('../outputs/string0.out')
#data = np.loadtxt('../data/minimum_2')
i=1
n=0
points = lx*lz
for n in range(0,N):
    fevalout = open("results/eigenvalues%dlz%dbylx%d.csv" % (n,lz,lx),'w')
    fevalout.write("x, y, z, eval\n")
    fevecout = open("results/eigenvectors%dlz%dbylx%d.csv" % (n,lz,lx),'w')
    fevecout.write("x, y, z, vx, vy, vz\n")
    fmod = open("results/modulusq%dlz%dbylx%d.csv" % (n,lz,lx),'w')
    fmod.write("x, y, z, vx, vy, vz\n")

    i=1
    print("start = ")
    print((i-1)*5 + 5*n*lx*lz)

    for i in range(0,(lz*lx)):

        indexz = np.mod(int(i), lz)
        indexx = int(int(i)/lz)

        Q1 = data[(i)*5 + 5*n*lx*lz]
        Q2 = data[((i)*5)+1 + 5*n*lx*lz]
        Q3 = data[((i)*5)+2 + 5*n*lx*lz]
        Q4 = data[((i)*5)+3 + 5*n*lx*lz]
        Q5 = data[((i)*5)+4 + 5*n*lx*lz]
        M = np.array([[Q1,Q2,Q3],[Q2,Q4,Q5],[Q3,Q5,-Q1-Q4]])

        mod = np.sqrt(Q1**2 + Q2**2 + Q3**2 + Q4**2 + Q5**2)
        #eval = np.array(np.linalg.eigvals(M))
        eval, evec = np.linalg.eig(M)
        idx = eval.argsort()[::-1]
        eval = eval[idx]
        evec = evec[:,idx]
        #print(i, eval, evec[:,])
        #i+=1

        fevalout.write("%d, 1, %d, %f\n" % (indexx, indexz, eval[0]))
        fmod.write("%d, 1, %d, %f\n" % (indexx, indexz, mod))
        fevecout.write("%d, 1, %d, %f, %f, %f\n" % (indexx, indexz, evec[0,0], evec[1,0],evec[2,0] ))

        #fevecout.write("%d, 1, %d, " + str(evec[1,0])+", "+str(evec[2,0])+", "+str(evec[0,0])+"\n" % (indexx, indexz))
        #if( i % 9999 == 0):

    print("end = ")
    print(((i)*5)+4 + 5*n*lx*lz)
    print(n)

    fmod.close()
    fevalout.close()
    fevecout.close()
#print(eigval)
