import os, string, sys
import subprocess
import math

from math import exp

R = 0.0019872156
# Ideal gas constant in kcal/mol*K
T = 298.15
# Temperature in K 

RT = R*T

i = 0
j = 0
k = 0
l = 0

os.chdir('BindEn')
#set i to the total number of unique receptor conformations
for i in range(10):
    if os.path.exists('binding_EnergyReceptorStemNameConformer_%d_final.txt' % i):
#Ki calculated in uM
	
        with open ('binding_EnergyReceptorStemNameConformer_%d_final.txt' % i) as infile, open ('inhibition_ConstantsReceptorStemNameConformerFinal_%d.txt' % i, 'w') as outfile:		
	   for line in infile:
		deltaG=float(line)
		#The following steps are to calculate the inhibition constant, Ki, in uM)
		#The equation RT = deltaG*lnKi is used, and rearranged to solve for Ki
		x = deltaG/RT
		Ki = exp(x)
		KiuM = Ki * 1000000
		outfile.write('%s\n' % KiuM)		
    else:
         print 'Path does not exist.'
         pass

			
print 'Done extracting inhibition constants'

