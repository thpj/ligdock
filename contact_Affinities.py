import os, string, sys
import math
import numpy
import numpy.linalg
import MDAnalysis
import MDAnalysis.analysis.contacts
import MDAnalysis.analysis.helanal
import scipy.io as sio

from math import sqrt

rootdir = '/PATH/TO/YOUR/dockedpdbs'

g = 0;
h = 0;

def computeAffinity(alist_i, alist_j):
    nij = 0;
    for atom_i in alist_i:
        for atom_j in alist_j:
            tmp_x = atom_i[0] - atom_j[0]
            tmp_y = atom_i[1] - atom_j[1]
            tmp_z = atom_i[2] - atom_j[2]
            d = math.sqrt(tmp_x**2 + tmp_y**2 + tmp_z**2)
            if d <= 4.0:
                nij = nij + 1;
    if nij > 0:
    	return nij;
    else:
    	return 0;

if __name__=='__main__':

#    k = 1;
#    for ts in u.trajectory[0::10]:
#        q = numpy.zeros((nRes, nRes)); 
#        fQ = open('../../Bcl2-Beclin/2ABO-apo/affinities/Aij_' + str(k) + '.txt', 'w');
#        for i in range(0, nRes):
#            for j in range(i+3, nRes):
#                selI = 'resid ' + str(i+1) + ' and not (name H*)';
#                selJ = 'resid ' + str(j+1) + ' and not (name H*)';
#                fi = u.selectAtoms(selI).coordinates(); #print len(fi);
#                fj = u.selectAtoms(selJ).coordinates(); #print len(fj);
#                q[i][j] = computeAffinity(fi, fj)/math.sqrt(len(fi)*len(fj));
#                q[j][i] = q[i][j];
#		print >> fQ, i, j, q[i][j] 
#        k = k+1;
#	fQ.close();
#        if (k%100 == 0):
#           print 'done ' + str(k) + '...' 
# """
        Q = numpy.zeros((52,20001));
        for subdir, dirs, files in os.walk(rootdir):
            g = g + 1;
            #set g to total number of directories
            if g <= 20000:
                #set h to total number of docking experiments performed
                for h in range(10):
                    os.chdir('DirectoryName%d' % g);
                    print 'Current directory is DirectoryName%d' % g;
                    if os.path.exists('/PATH/TO/YOUR/DirectoryName%d/FinalPDBFileName%d_%d.pdb' % (g, g, h)):
                            print 'Path FinalPDBFileName exists!';
                            nRes = 52;
                            k = 1;
                            u = MDAnalysis.Universe('FinalPDBFileName%d_%d.pdb' % (g, h), permissive=True);
                            for ts in u.trajectory:
                                for i in range(0, nRes):
                                    sel53 = 'resid ' + str(53) + ' and not (name H*)';
                                    selJ = 'resid ' + str(i+1) + ' and backbone and (name N*)';
                                    f53 = u.selectAtoms(sel53).coordinates(); #print len(fi); 
                                    #print 'f53', f53.atoms
                                    fj = u.selectAtoms(selJ).coordinates(); #print len(fj);
                                    #print 'fj', fj.atoms
                                    Q[i, '%d' % g] = ((computeAffinity(f53, fj))/(sqrt(len(f53) * len(fj))));
                                k = k+1;
                                if k%100 == 0:
                                    print 'done ' + str(k) + '...';
                            os.chdir('../');
                    else:
		      print 'Path FinalPDBFileName does not exist!';
		      os.chdir('../');
		      continue;
            elif g > 20000:
	      print 'All done';
	      break;
        sio.savemat("/PATH/TO/SAVE/ContactsBackboneNligand.mat", { 'Q': Q });
#	q = q/2000.;
#	
#	fOut = open('AvgQ_P27.txt', 'w');
#	for i in range (0,52):
#		for j in range(0,52):
#			print >> fOut, q[i][j],
#		print >> fOut, '\n';
#	fOut.close();	
#
#	print q;

