import os, string, sys
import subprocess
import itertools
import csv

rootdir = '/PATH/TO/YOUR/dockingfiles'

i = 0 
j = 0 
h = 0
k = 0
l = 0
m = 0
n = 0
p = 0
q = 0

x = 0 
y = 0 
w = 0 
v = 0 

# Save ligand and receptor for each docking experiment as a PDBQT file, even if docking is not the lowest energy.
# write_all_complexes.py is used to do this. Although _(some number) is added to the end of the file name,
# this does not indicate cluster rank. It only indicates this was the conformer from search 0, 1, 2, etc...
for subdir, dirs, files in os.walk(rootdir):
    i = i + 1
# Set the following condition to total number of ligands or directories where docking log files are located
    if i <= 63:
        os.chdir('Z_%d' % i)
        print 'Z_%d_toPDBQT' % i
        for j in range(10):
            if os.path.exists('YourDockingLog%d_%d.dlg' % (j, i)):
                os.system('/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/write_all_complexes.py -d YourDockingLogFile%d_%d.dlg -r YourReceptorFile%d.pdbqt -o PATH/TO/YOUR/dockedpdbs/PDBFileStemName%d_%d' % (j, i, j, j, i))
            else:
                print 'Docking log does not exist.'
                pass
        os.chdir('../')
    elif i > 63:
        break 

# Remove the Q and T columns from the PDBQT file  
os.chdir('/PATH/TO/YOUR/dockedpdbs')
# Set x range to number of ligands or directories from docking experiment
for x in range(1,64):
    # Set y number of searches performed in docking experiment. For example, we generate 10 runs, so only 10 docking experiments.
    # The range is determined from the previously used write_all_complexes.py script, which starts at 0 until the last pdbqt file 
    # is written.
    for y in range(10):
        # Set w range to total number of unique receptor conformers used in the docking experiments, based off of some base_filename
        # structure. In our case, we had no more than 10 receptor conformers linked to a unique base_filename.
        for w in range(10):
            if os.path.exists('PDBFileStemName%d_%d_%d.pdbqt' % (w, x, y)):
                print 'Path does exist.'
                convert_PDBQTtoPDB = 'cut -c-66 PDBFileStemName%d_%d_%d.pdbqt > PDBFileStemName%d_%d_%d.pdb' % (w, x, y, w, x, y)
                subprocess.call(convert_PDBQTtoPDB, shell=True)
                os.remove('PDBFileStemName%d_%d_%d.pdbqt' % (w, x, y))
                v = v + 1
                print 'v', v
                # Change x, y, and w in order to reflect your file structure in order to break the loop, which removes the q and t
                # columns from the pdbqt to make it a pdb file. 
                if x < 8 and y < 9 and w < 9:
                    continue
                elif x == 8 and y == 9 and w == 9:
                    break
            else:
                print 'Path does not exist.'
                pass

h = 0
k = 0
l = 0

q = 0
r = 0
s = 0
t = 0
u = 0
v = 0
w = 0
x = 0
a = 0
# Set this range to number of ligands or directories from docking experiment
for h in range(1,30440):
    # set to number of searches performed in docking experiment. For example, we generate 10 runs, so only 10 docking experiments.
    # the range is determined from the previously used write_all_complexes.py script, which starts at 0 until the last pdbqt file 
    # is written.
    for k in range(10):
        # set range to total number of unique receptor conformers used in the docking experiments, based off of some base_filename
        #structure. In our case, we had no more than 10 receptor conformers linked to a unique base_filename.
        for l in range(0,1):
            if os.path.exists('PDBFileStemName%d_%d_%d.pdb' % (l, h, k)):
                print 'PDB file exists, number %d.' % h
                fout1 = open('PDBFileStemName%d_%d_%d_final.pdb' % (l, h, k), 'w')
                with open('PDBFileStemName%d_%d_%d.pdb' % (l, h, k), 'r') as fin1:
                    # The following set of conditions renames that atom types of the ligand only, such as C to C1, by each instance it appears.
                    # It also renumbers the integer to be consistent with the number of residues. The second string insertion while writing the line
                    # will need to be changed according to the protein under study. Current implementation only renumbers atom types of C, N, O, H
                    # F, S, P, Cl, and Br. ***Important*** Only formatted for C, N, O, H, F, S, and P to have more than 999 atoms; Br and Cl only for 99 atoms.
                    # If the following hard conditions are not met, then the line is just written as is.
                    for line in fin1:
                        if line[17:20]=='LIG' and line [13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:13] + 'C%d' % q + line[15:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:13] + 'C%d' % q + line[16:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + line[16:23] + '314 ' + line[27:80])
                        elif line[17:20]=='LIG' and line [13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:13] + 'N%d' % r + line[15:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:13] + 'N%d' % r + line[16:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + line[16:23] + '314 ' + line[27:80])
                        elif line[17:20]=='LIG' and line [13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:13] + 'O%d' % s + line[15:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:13] + 'O%d' % s + line[16:23] + '314 ' + line[27:80])
                            elif s >=100:
                                fout1.write(line[0:12] + 'O%d' % s + line[16:23] + '314 ' + line[27:80])
                        elif line[17:20]=='LIG' and line [13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:13] + 'H%d' % t + line[15:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:13] + 'H%d' % t + line[16:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + line[16:23] + '314 ' + line[27:80])
                        elif line[17:20]=='LIG' and line [13:14] == 'F':
                            u = u + 1
                            if u < 10:
                                fout1.write(line[0:13] + 'F%d' % u + line[15:23] + '314 ' + line[27:80])
                            elif u < 100:
                                fout1.write(line[0:13] + 'F%d' % u + line[16:23] + '314 ' + line[27:80])
                            elif u >= 100:
                                fout1.write(line[0:12] + 'F%d' % u + line[16:23] + '314 ' + line[27:80])
                        elif line[17:20]=='LIG' and line [13:14] == 'S':
                            v = v + 1
                            if v < 10:
                                fout1.write(line[0:13] + 'S%d' % v + line[15:23] + '314 ' + line[27:80])
                            elif v < 100:
                                fout1.write(line[0:13] + 'S%d' % v + line[16:23] + '314 ' + line[27:80])
                            elif v >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + line[16:23] + '314 ' + line[27:80])
                        elif line[17:20]=='LIG' and line [12:14] == 'BR': 
                            w = w + 1
                            if w < 10:
                                fout1.write(line[0:12] + 'BR%d' % w + line[16:23] + '314 ' + line[27:80])
                            elif w >= 10:
                                fout1.write(line[0:12] + 'BR%d' % w + line[16:23] + '314' + line[27:80])
                        elif line[17:20]=='LIG' and line [12:14] == 'CL':
                            x = x + 1
                            if x < 10:
                                fout1.write(line[0:12] + 'CL%d' % x + line[16:23] + '314 ' + line[27:80])
                            elif x >= 10:
                                fout1.write(line[0:12] + 'CL%d' % x + line[16:23] + '314' + line[27:80])
                        elif line[17:20]=='LIG' and line [13:14] == 'P':
                            a = a + 1
                            if a < 10:
                                fout1.write(line[0:13] + 'P%d' % a + line[15:23] + '314 ' + line[27:80])
                            elif a < 100:
                                fout1.write(line[0:13] + 'P%d' % a + line[16:23] + '314 ' + line[27:80])
                            elif a >= 100:
                                fout1.write(line[0:12] + 'P%d' % a + line[16:23] + '314 ' + line[27:80])
                        elif line[17:20]=='UNK' and line [13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:13] + 'C%d' % q + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:13] + 'C%d' % q + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >= 100:
                                fout1.write(line[0:12] + 'C%d' % q + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='UNK' and line [13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:13] + 'N%d' % r + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:13] + 'N%d' % r + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='UNK' and line [13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:13] + 'O%d' % s + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:13] + 'O%d' % s + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='UNK' and line [13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:13] + 'H%d' % t + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:13] + 'H%d' % t + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='UNK' and line [13:14] == 'F':
                            u = u + 1
                            if u < 10:
                                fout1.write(line[0:13] + 'F%d' % u + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif u < 100:
                                fout1.write(line[0:13] + 'F%d' % u + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif u >= 100:
                                fout1.write(line[0:12] + 'F%d' % u + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='UNK' and line [13:14] == 'S':
                            v = v + 1
                            if v < 10:
                                fout1.write(line[0:13] + 'S%d' % v + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif v < 100:
                                fout1.write(line[0:13] + 'S%d' % v + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif v >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='UNK' and line [12:14] == 'BR': 
                            w = w + 1
                            if w < 10:
                                fout1.write(line[0:12] + 'BR%d' % w + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif w >= 10:
                                fout1.write(line[0:12] + 'BR%d' % w + line[16:17] + 'LIG' + line[20:23] + '314' + line[27:80])
                        elif line[17:20]=='UNK' and line [12:14] == 'CL':
                            x = x + 1
                            if x < 10:
                                fout1.write(line[0:12] + 'CL%d' % x + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif x >= 10:
                                fout1.write(line[0:12] + 'CL%d' % x + line[16:17] + 'LIG' + line[20:23] + '314' + line[27:80])
                        elif line[17:20]=='UNK' and line [13:14] == 'P':
                            a = a + 1
                            if a < 10:
                                fout1.write(line[0:13] + 'P%d' % a + line[15:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif a < 100:
                                fout1.write(line[0:13] + 'P%d' % a + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif a >= 100:
                                fout1.write(line[0:12] + 'P%d' % a + line[16:17] + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='<1>' and line [13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:13] + 'C%d' % q + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:13] + 'C%d' % q + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >= 100:
                                fout1.write(line[0:12] + 'C%d' % q + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='<1>' and line [13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:13] + 'N%d' % r + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:13] + 'N%d' % r + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='<1>' and line [13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:13] + 'O%d' % s + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:13] + 'O%d' % s + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='<1>' and line [13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:13] + 'H%d' % t + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 10:
                                fout1.write(line[0:13] + 'H%d' % t + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >=100:
                                fout1.write(line[0:12] + 'H%d' % t + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='<1>' and line [13:14] == 'F':
                            u = u + 1
                            if u < 10:
                                fout1.write(line[0:13] + 'F%d' % u + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif u < 100:
                                fout1.write(line[0:13] + 'F%d' % u + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif u >= 100:
                                fout1.write(line[0:12] + 'F%d' % u + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='<1>' and line [13:14] == 'S':
                            v = v + 1
                            if v < 10:
                                fout1.write(line[0:13] + 'S%d' % v + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif v < 100:
                                fout1.write(line[0:13] + 'S%d' % v + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif v >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:20]=='<1>' and line [12:14] == 'BR': 
                            w = w + 1
                            if w < 10:
                                fout1.write(line[0:12] + 'BR%d' % w + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif w >= 10:
                                fout1.write(line[0:12] + 'BR%d' % w + ' ' + 'LIG' + line[20:23] + '314' + line[27:80])
                        elif line[17:20]=='<1>' and line [12:14] == 'CL':
                            x = x + 1
                            if x < 10:
                                fout1.write(line[0:12] + 'CL%d' % x + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif x >= 10:
                                fout1.write(line[0:12] + 'CL%d' % x + ' ' + 'LIG' + line[20:23] + '314' + line[27:80])
                        elif line[17:20]=='<1>' and line [13:14] == 'P':
                            a = a + 1
                            if a < 10:
                                fout1.write(line[0:13] + 'P%d' % a + '  ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif a < 100:
                                fout1.write(line[0:13] + 'P%d' % a + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                            elif a >= 100:
                                fout1.write(line[0:12] + 'P%d' % a + ' ' + 'LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLY d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLY d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLY d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLY d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLY d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ALA d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ALA d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ALA d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ALA d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ALA d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'VAL d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'VAL d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'VAL d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'VAL d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'VAL d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LEU d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LEU d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LEU d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LEU d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LEU d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ILE d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ILE d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ILE d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ILE d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ILE d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'MET d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'MET d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'MET d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'MET d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'MET d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PRO d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PRO d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PRO d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PRO d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PRO d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PHE d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PHE d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PHE d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PHE d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'PHE d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TRP d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TRP d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TRP d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TRP d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TRP d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'SER d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'SER d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'SER d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'SER d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'SER d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'THR d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'THR d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'THR d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'THR d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'THR d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASN d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASN d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASN d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASN d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASN d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLN d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLN d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLN d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLN d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLN d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TYR d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TYR d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TYR d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TYR d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'TYR d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'CYS d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'CYS d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'CYS d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'CYS d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'CYS d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LYS d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LYS d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LYS d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LYS d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'LYS d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ARG d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ARG d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ARG d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ARG d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ARG d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'HIS d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'HIS d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'HIS d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'HIS d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'HIS d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASP d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASP d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASP d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASP d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'ASP d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLU d' and line[13:14] == 'C':
                            q = q + 1
                            if q < 10:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q < 100:
                                fout1.write(line[0:12] + ' ' + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif q >=100:
                                fout1.write(line[0:12] + 'C%d' % q + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLU d' and line[13:14] == 'N':
                            r = r + 1
                            if r < 10:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r < 100:
                                fout1.write(line[0:12] + ' ' + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif r >= 100:
                                fout1.write(line[0:12] + 'N%d' % r + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLU d' and line[13:14] == 'O':
                            s = s + 1
                            if s < 10:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s < 100:
                                fout1.write(line[0:12] + ' ' + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif s >= 100:
                                fout1.write(line[0:12] + 'O%d' % s + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLU d' and line[13:14] == 'H':
                            t = t + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'H%d' % t + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        elif line[17:22] == 'GLU d' and line[13:14] == 'S':
                            v = v + 1
                            if t < 10:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + '  LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t < 100:
                                fout1.write(line[0:12] + ' ' + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                            elif t >= 100:
                                fout1.write(line[0:12] + 'S%d' % v + ' LIG' + line[20:23] + '314 ' + line[27:80])
                        else:
                            fout1.write(line)
                fin1.close()
                fout1.close()
                # Removes initially created pdb
                #os.remove('PDBFileStemName%d_%d_%d.pdb' % (l, h, k))
                q = 0
                r = 0
                s = 0
                t = 0
                u = 0
                v = 0
                w = 0
                a = 0
                x = 0
                if h < 8 and k < 9 and l < 9:
                    continue
                elif h == 8 and k == 9 and l == 9:
                    break
            else:
                print 'PDB file does not exist!'
                continue

h = 0
k = 0
l = 0
#set z to number of atoms in protein/receptor
z = 3066 
# Set this range to number of ligands or directories from docking experiment
for h in range(1,30440):
    # Set to number of searches performed in docking experiment. For example, we generate 10 runs, so only 10 docking experiments.
    # The range is determined from the previously used write_all_complexes.py script, which starts at 0 until the last pdbqt file 
    # is written.
    for k in range(10):
        # Set range to total number of unique receptor conformers used in the docking experiments, based off of some base_filename
        # structure. In our case, we had no more than 10 receptor conformers linked to a unique base_filename.
        for l in range(0,1):
            if os.path.exists('PDBFileStemName%d_%d_%d_final.pdb' % (l, h, k)):
                print 'Final PDB exists, %d!' % h
                fout2 = open('PDBFileStemName%d_%d_%d_final2.pdb' % (l, h, k), 'w')
                with open('PDBFileStemName%d_%d_%d_final.pdb' % (l, h, k), 'r') as fin2:
                    # Changes the atom serial number to be consistent with the receptor protein atom serial numbers
                    # In this instance, the final number 3066. However, this will need to be changed based on the
                    # the protein used. This is controlled by the variable 'z'.
                    for line in fin2:
                        if line[17:20] == 'LIG' and line[0:5] == 'ATOM ':
                            z = z + 1
                            fout2.write(line[0:7] + '%d ' % z + line[12:80])
                    # For some reason or another, the TER line was not included when the PDBQT file was written...possibly not in original file.
                    # The following condition adds the appropriate TER line to indicate separation between the ligand and protein.
                    #    elif line [0:26] == 'ATOM   3065  OXT PHE A 313':
                    #        fout2.write(line + 'TER    3066      PHE A 313                                        ''\n')
                        else:
                            fout2.write(line)
                fin2.close()
                fout2.close()
                # Removes the first created *_final.pdb
                #os.remove('PDBFileStemName%d_%d_%d_final.pdb' % (l, h, k))
                #***Don't forget to reset z here!!!***
                z = 3066
                if h < 8 and k < 9 and l < 9:
                    continue
                elif h == 8 and k == 9 and l == 9:
                    break
                # Final PDB will have the following naming structure. This is only given as an example.
                # Example filename PDB9_8_9_final2.pdb
print 'Final PDB created'