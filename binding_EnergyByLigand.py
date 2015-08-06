import os, string, sys
import subprocess

rootdir = 'PATH/TO/YOUR/completedDocking'

i = 0
j = 0
h = 0
k = 0
l = 0
m = 0
n = 0
p = 0
q = 0
r = 0

#Extracts a text file with clustering information, RMSD, binding energy, and other information summarizing the docking experiment by using summarize_docking.py.
#-t flag is set to zero to make RMSD tolerance zero. Otherwisre, summarize_docking.py only writes out the information pertaining to the lowest energy conformer
#in each cluster.
#-a flag is set to append the next set of docking summaries to the end of the previous file.
#The script is divided into four to obtain this information based on some base filename structure.

for subdir, dirs, files in os.walk(rootdir):
    i = i + 1
# set i for total number of ligands
    if i <= 200:
        os.chdir('Z_%d' % i)
        print 'Z_%d' % i
        #set j to the total number of unique receptor conformations
        for j in range(10):
            if os.path.exists('dockingLogStemName%d_%d.dlg' % (j, i)):
                print 'Docking log exists'
                os.system('/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/summarize_docking.py -l dockingLogStemName%d_%d.dlg -t 0 -a -o /PATH/TO/WRITE/YOUR/BindEn/binding_Energyligand_%d.txt' % (j, i, i))
            else:
                print 'Docking log does not exist'
                continue
        os.chdir('../')
    elif i > 200:
        print 'Done summarizing docking of dockingLogStemName%d_%d' % (j, i)
        break
                
#Since we are only interested in obtaining the binding energy for this part, the next step extracts the binding energy column only.
       
os.chdir('BindEn')
for q in range(201):
#set q for total number of ligands
         if os.path.exists('binding_Energyligand_%d.txt' % q):
             fout = open('binding_Energyligand_%d_2.txt' % q, 'a')
             with open('binding_Energyligand_%d.txt' % q) as fin:
                 data=fin.read().splitlines(True)
                 fout.writelines(data[1:])
                 #This part removes the header information created in the text file generated from summarize_docking.py
         else:
             print 'Path does not exist.'
             continue
#set r to number of ligands
for r in range(201):
         if os.path.exists('binding_Energyligand_%d_2.txt' % r):
             extract_EnergyLigand = 'cut -d , -f 3 binding_Energyligand_%d_2.txt > binding_Energyligand_%d_final.txt' % (r, r)
             subprocess.call(extract_EnergyLigand, shell =True)
             #This part extracts only the third column, using ',' delimiters to grab the third column.
         else:
             print 'Path does not exist'
             continue

print 'Done extracting binding energies'