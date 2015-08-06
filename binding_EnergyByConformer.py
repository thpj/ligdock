import os, string, sys
import subprocess

rootdir = '/PATH/TO/YOUR/completedDocking'

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
s = 0
t = 0

#Extracts a text file with clustering information, RMSD, binding energy, and other information summarizing the docking experiment by using summarize_docking.py.
#-t flag is set to zero to make RMSD tolerance zero. Otherwisre, summarize_docking.py only writes out the information pertaining to the lowest energy conformer
#in each cluster.
#-a flag is set to append the next set of docking summaries to the end of the previous file.

for subdir, dirs, files in os.walk(rootdir):
#set i to the total number of ligands or directories
    i = i + 1
    if i <= 200:
        os.chdir('Z_%d' % i)
        print 'Z_%d' % i
        #set j to total number of unique receptor conformations
        for j in range(10):
            if os.path.exists('dockingLogStemName%d_%d.dlg' % (j, i)):
                print 'Docking log exists'
                os.system('/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/summarize_docking.py -l dockingLogStemName%d_%d.dlg -t 0 -a -o /PATH/TO/WRITE/YOUR/BindEn/binding_EnergyReceptorStemNameConformer_%d.txt' % (j, i, j))
            else:
                print 'Docking log does not exist'
                continue
        os.chdir('../')
    elif i > 200:
        print 'Done summarizing docking of dockingLogStemName%d_%d' % (j, i)
        break
        
#Since we are only interested in obtaining the binding energy for this part, the next step extracts the binding energy column only. 
       

os.chdir('BindEn')
#set for total number of uniquely identified receptor conformations
for q in range(10): 
         if os.path.exists('binding_EnergyReceptorStemNameConformer_%d.txt' % q):
            fout = open('binding_EnergyReceptorStemNameConformer_%d_2.txt' % q, 'a')
            with open('binding_EnergyReceptorStemNameConformer_%d.txt' % q) as fin:
                data= fin.read().splitlines(True)
                fout.writelines(data[1:])
#this step removes any of the headers from the initial file, such as 'BE', 'Rank','Cluster', etc... generated from summarize_docking.py
         else:
            print 'Path does not exist.'
            continue
#set r to number of unique receptor conformations          
for r in range(10):
            if os.path.exists('binding_EnergyReceptorStemNameConformer_%d_2.txt' % r):
                extract_EnergyLigand01 = 'cut -d , -f 3 binding_EnergyReceptorStemNameConformer_%d_2.txt > binding_EnergyReceptorStemNameConformer_%d_final.txt' % (r, r)
                subprocess.call(extract_EnergyLigand01, shell=True)
#this step extracts only the third column, the binding energies, from the text file without commas. No further cleanning up of the data 
#should be needed at this point.
            else:
                print 'Path does not exist.'
                continue     
print 'Done extracting binding energies'