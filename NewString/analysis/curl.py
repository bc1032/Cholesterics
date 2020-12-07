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

    fcurl = open("curl/curlq%dlz%dbylx%d.csv" % (n,lz,lx),'w')
    fcurl.write("#x, y, z, vx, vy, vz\n")

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
        curlq = 0.0
        M = np.array([[Q1,Q2,Q3],[Q2,Q4,Q5],[Q3,Q5,-Q1-Q4]])

        if i <= 5*lz*lz - 5*lz and i >= 5*lz:
            if (i)*5 + 5*lx*lz % 5*lz == 0:
                Q1L = data[(i)*5 + 5*n*lx*lz + 5*lx - 5]
                Q2L = data[((i)*5)+1 + 5*n*lx*lz + 5*lx - 5]
                Q3L = data[((i)*5)+2 + 5*n*lx*lz + 5*lx - 5]
                Q4L = data[((i)*5)+3 + 5*n*lx*lz + 5*lx - 5]
                Q5L = data[((i)*5)+4 + 5*n*lx*lz + 5*lx - 5]
            elif ((i)*5 + 5*lx*lz) % 5*lz == 1:
                Q1R = data[(i)*5 + 5*n*lx*lz - 5*lz + 5]
                Q2R = data[((i)*5)+1 + 5*n*lx*lz - 5*lx + 5]
                Q3R = data[((i)*5)+2 + 5*n*lx*lz - 5*lx + 5]
                Q4R = data[((i)*5)+3 + 5*n*lx*lz - 5*lx + 5]
                Q5R = data[((i)*5)+4 + 5*n*lx*lz - 5*lx + 5]
            else:
                Q1L = data[(i)*5 + 5*n*lx*lz - 5]
                Q2L = data[((i)*5)+1 + 5*n*lx*lz - 5]
                Q3L = data[((i)*5)+2 + 5*n*lx*lz - 5]
                Q4L = data[((i)*5)+3 + 5*n*lx*lz - 5]
                Q5L = data[((i)*5)+4 + 5*n*lx*lz - 5]

                Q1R = data[(i)*5 + 5*n*lx*lz + 5]
                Q2R = data[((i)*5)+1 + 5*n*lx*lz + 5]
                Q3R = data[((i)*5)+2 + 5*n*lx*lz + 5]
                Q4R = data[((i)*5)+3 + 5*n*lx*lz + 5]
                Q5R = data[((i)*5)+4 + 5*n*lx*lz + 5]

            Q1A = data[(i)*5 + 5*n*lx*lz + 5*lz]
            Q2A = data[((i)*5)+1 + 5*n*lx*lz + 5*lz]
            Q3A = data[((i)*5)+2 + 5*n*lx*lz + 5*lz]
            Q4A = data[((i)*5)+3 + 5*n*lx*lz + 5*lz]
            Q5A = data[((i)*5)+4 + 5*n*lx*lz + 5*lz]

            Q1B = data[(i)*5 + 5*n*lx*lz - 5*lz]
            Q2B = data[((i)*5)+1 + 5*n*lx*lz - 5*lz]
            Q3B = data[((i)*5)+2 + 5*n*lx*lz - 5*lz]
            Q4B = data[((i)*5)+3 + 5*n*lx*lz - 5*lz]
            Q5B = data[((i)*5)+4 + 5*n*lx*lz - 5*lz]

            MA = np.array([[Q1A,Q2A,Q3A],[Q2A,Q4A,Q5A],[Q3A,Q5A,-Q1A-Q4A]])
            MB = np.array([[Q1B,Q2B,Q3B],[Q2B,Q4B,Q5B],[Q3B,Q5B,-Q1B-Q4B]])
            ML = np.array([[Q1L,Q2L,Q3L],[Q2L,Q4L,Q5L],[Q3L,Q5L,-Q1L-Q4L]])
            MR = np.array([[Q1R,Q2R,Q3R],[Q2R,Q4R,Q5R],[Q3R,Q5R,-Q1R-Q4R]])

            lciv = np.array([[0,1,-1],[1,0,1],[-1,1,0]])
            MZ = MA - MB
            MX = MR - ML

            curlq = 0.0
            for qi in range(0,3):
                for qj in range(0,3):
                    for ql in range(0,3):
                        curlq = curlq + M[qi][qj]*MZ[ql][qi]*lciv[qj][ql]
                        curlq = curlq + M[qi][qj]*MX[ql][qi]*lciv[qj][ql]
            #print(curlq)
        fcurl.write("%f\n" % (curlq))

    print("end = ")
    print(((i)*5)+4 + 5*n*lx*lz)
    print(n)

    fcurl.close()
