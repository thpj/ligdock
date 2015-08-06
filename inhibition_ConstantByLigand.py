import os, string, sys
import subprocess
import math

from math import exp

R = 0.0019872156
# Ideal gas constant in kcal/mol*K
T = 298.15
# Temperature in K 

RT = R*T 
#in kcal/mol
i = 0
os.chdir('BindEn')
#change directory to where your binding energies are
#set i to the number of ligands
for i in range(201):
    if os.path.exists('binding_Energyligand_%d_final.txt' % i):
#Ki calculated in uM
		
        with open ('binding_Energyligand_%d_final.txt' % i) as infile, open ('inhibition_ConstantsLigand_%d_final.txt' % i, 'w') as outfile:		
	   for line in infile:
		deltaG=float(line)
		#The following steps are to calculate the inhibition constant, Ki, in uM)
		#The equation RT = deltaG*lnKi is used, and rearranged to solve for Ki
		x = deltaG/R
		Ki = exp(x)
		KiuM = Ki * 1000000
		outfile.write('%s\n' % KiuM)
					
print 'Done extracting inhibition constants'