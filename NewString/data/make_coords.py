import numpy as np
import pickle

GRIDX = 30
GRIDY = 1
GRIDZ = 30
w_1 = 4
w_2 = 2
WIDTH=GRIDZ-1
PI = np.arctan(1)*4.0
theta = 0.0
Q1, Q2, Q3, Q4, Q5 = 0.0, 0.0, 0.0, 0.0, 0.0
PHIS = np.empty((GRIDX,GRIDY,GRIDZ))

A = 0.3
B = 0.2
C = 1.0

s = (B + np.sqrt(B**2 + 24*A*C))/(4.0*C)
#long_array(:)=0.0
#PHIS(:,:,:) = -1.0

with open('minimum_1', 'w') as output_file:
	for J1 in range(0,GRIDY):
		for J2 in range(0,GRIDZ):
			for J3 in range(0,GRIDX):
				theta = w_1*PI*(J3)/WIDTH
				Q1 = s*(np.cos(theta)*np.cos(theta) - (1.0/3.0))
				Q2 = s*(np.cos(theta)*np.cos(theta))
				Q3 = 0.0
				Q4 = s*(np.sin(theta)*np.sin(theta) - (1.0/3.0))
				Q5 = 0.0
				output_file.write(str(Q1)+"\n")
				output_file.write(str(Q2)+"\n")
				output_file.write(str(Q3)+"\n")
				output_file.write(str(Q4)+"\n")
				output_file.write(str(Q5)+"\n")
with open('minimum_2', 'w') as output_file:
	for J1 in range(0,GRIDY):
		for J2 in range(0,GRIDZ):
			for J3 in range(0,GRIDX):
				theta = w_2*PI*(J3)/WIDTH
				Q1 = s*(np.cos(theta)*np.cos(theta) - (1.0/3.0))
				Q2 = s*(np.cos(theta)*np.cos(theta))
				Q3 = 0.0
				Q4 = s*(np.sin(theta)*np.sin(theta) - (1.0/3.0))
				Q5 = 0.0
				output_file.write(str(Q1)+"\n")
				output_file.write(str(Q2)+"\n")
				output_file.write(str(Q3)+"\n")
				output_file.write(str(Q4)+"\n")
				output_file.write(str(Q5)+"\n")
