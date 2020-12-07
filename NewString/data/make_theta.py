import numpy as np
import pickle

GRIDX = 10
GRIDY = 1
GRIDZ = 10
w_1 = 4
w_2 = 2
WIDTH=GRIDZ-1
PI = np.arctan(1)*4.0


PHIS = np.empty((GRIDX,GRIDY,GRIDZ))
#long_array(:)=0.0
#PHIS(:,:,:) = -1.0


for J1 in range(0,GRIDX):
	for J2 in range(0,GRIDY):
		for J3 in range(0,GRIDZ):
			PHIS[J1,J2,J3] = (w_2*PI*(J3)/WIDTH)
print(PHIS)

with open('theta_1', 'w') as output_file:
	for J1 in range(0,GRIDY):
		for J2 in range(0,GRIDZ):
			for J3 in range(0,GRIDX):
				output_file.write( str(w_1*PI*(J3)/WIDTH)+"\n")

with open('theta_2', 'w') as output_file:
	for J1 in range(0,GRIDY):
		for J2 in range(0,GRIDZ):
			for J3 in range(0,GRIDX):
					output_file.write( str(w_2*PI*(J3)/WIDTH)+"\n")
