import os, string, sys
import subprocess
import MDAnalysis
import numpy

from MDAnalysis import *



rootdir = '/PATH/TO/YOUR/dockedpdbs';

e = 0;
b = 0;
t = 0;
h = 0;
p = 0;
j = 0;

u = 0;
w = 0
os.chdir('PATH/TO/YOUR/dockedpdbs/')
#set e range to total number of ligands used in your docking experiment
for e in range(1,30440):
    # set u to the total number of unique receptor conformations/IDs
	for u in range(1):
	    # set w to the total number of docking experiments carried out with each ligand to receptor conformation. By default,
	    #AutoDock4 performs 10 dockings in each experiment, so this number will usually be 10, unless your specify otherwise.
            for w in range(10):
                if os.path.exists('PDBFileStemName%d_%d_%d_final2.pdb' % (u, e, w)):
                    #create cent_Masspdbs directory before starting
                    outfile = open('cent_Masspdbs/PDBFileStemName%d_%d_CentMass.pdb' % (u, e,), 'a')
                    print 'Path to final PDB exists!'
                    buf = ''
                    #load your docked PDB into the MDAnalysis universe
                    u1 = MDAnalysis.Universe('PDBFileStemName%d_%d_%d_final2.pdb' % (u, e, w), permissive=False)
                    print 'Got the PDBFileStemName%d_%d_%d_final2.pdb' % (u, e, w)
                    #set resid to the residue number of interest of your ligand in the docked PDB file.
                    LIG = u1.selectAtoms('resid 3')
                    print 'LIG', LIG
                    CentMass =LIG.centerOfMass()
                    #CentMass calculates center of mass
                    print 'Center of mass is:', CentMass;
                    #CentMassStr returns the CentMass array as a string
                    CentMassStr = numpy.array_str(CentMass)
                    StripCM = CentMassStr.strip('[]')
                    #StripCM removes the brackets so that the string can be written to the final center of mass PDB.
                    print 'String center of mass is:', StripCM
                    with open('PDBFileStemName%d_%d_%d_final2.pdb' % (u, e, w)) as infile:
                        for line in infile:
                            #set int(line[7:11] equivalent to the first atom number of your ligand and line[17:20] to the residue ID LIG or something else if you have chosen a different ID.
                            if int(line[7:11]) == 1500 and line[17:20] == 'LIG':
                                buf = line[0:5]+ '  %d' % (1500 +  w)+ '  CA ' + line[16:23] + '%d' % (3 + w) + '      '+ StripCM[1:7] + ' '+StripCM[13:20]+ '  '+StripCM[27:33]+  '  0.00' + '  0.00' '\n' 'TER' '\n';
                                outfile.write(buf)
                                #setting the below condition only writes the protein structure once to the file by setting the condition of w to zero. Otherwise, the protein will be written
                                #for each docking experiment done with the same receptor/ligand combination into one file.
			    elif int(line[7:11]) < 1500 and w == 0:
				outfile.write(line)
                            else:
                                continue
                        infile.close()
                    outfile.close()
                else:
                    print 'Final PDB file does not exist'
                    continue
print 'Finished files'
