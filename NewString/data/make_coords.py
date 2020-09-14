import numpy as np
import pickle

GRIDX = 1
GRIDY = 30
GRIDZ = 30
w = 4
WIDTH=GRIDZ-1
PI = np.arctan(1)*4.0


PHIS = np.empty((GRIDX,GRIDY,GRIDZ))
#long_array(:)=0.0
#PHIS(:,:,:) = -1.0


for J1 in range(0,GRIDX):
	for J2 in range(0,GRIDY):
		for J3 in range(0,GRIDZ):
			PHIS[J1,J2,J3] = (w*PI*(J3)/WIDTH)
print(PHIS)

with open('minimum_1', 'w') as output_file:
	for J1 in range(0,GRIDX):
		for J2 in range(0,GRIDY):
			for J3 in range(0,GRIDZ):
				if J3 < GRIDZ-1:
					#PHIS[J1,J2,J3] = (w*PI*(J3)/WIDTH)
					output_file.write( str(w*PI*(J3)/WIDTH)+"\n")
					#output_file.write( str(w*PI*(J3)/WIDTH) + "	" + str(w*PI*(J3)/WIDTH) + "	" + str(0)  +"\n")
				else:
					#PHIS[J1,J2,J3] = 1.0
					output_file.write( str(w*PI*(J3)/WIDTH)+"\n")
					#output_file.write( str(w*PI) + "	" + str(w*PI) + "	" + str(0.0)  +"\n")
				#output_file.write( str(w*PI*(J3)/WIDTH) + "	" + str(w*PI*(J3)/WIDTH) + "	" + str(w*PI*(J3)/WIDTH)  +"\n")
#data = np.loadtxt('coordspy.dat')
#with open('coords', 'w') as newfile:
#	for line in range(0,len(data)):
#		newfile.write( str(data[line,0]) + "\n" + str(data[line,1]) + "\n" + str(data[line,2])  +"\n")
			
#with open('cossincoordspy2.dat', 'w') as output_file:
#	for J1 in range(0,GRIDX):
#		for J2 in range(0,GRIDY):
#			for J3 in range(0,GRIDZ):
#				if J3 < GRIDZ-1:
				#PHIS[J1,J2,J3] = (w*PI*(J3)/WIDTH)
#					output_file.write( str(np.cos(w*PI*(J3)/WIDTH)) + "	" + str(np.sin(w*PI*(J3)/WIDTH)) + "	" + str(0)+ "	" + str(J3)+ "	" + str(J2)+ "	" + str(J1)  +"\n")
#				else:
#					output_file.write( str(1.0) + "	" + str(0) + "	" + str(0)+ "	" + str(J3)+ "	" + str(J2)+ "	" + str(J1)  +"\n")
				#output_file.write( str(w*PI*(J3)/WIDTH) + "	" + str(w*PI*(J3)/WIDTH) + "	" + str(w*PI*(J3)/WIDTH)  +"\n")
		
