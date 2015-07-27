import os, string, sys, glob
import subprocess

i = 0
#Unless your computer has losts of disk space, you're going to want to remove the
#map, fld, and xyz files used in generating your grids and after your docking
#experiment is complete.

#adjust i such that it reflect the directories you wish to remove files from
for i in range(64, 101):
    if os.path.exists('Z_%d' % i):
        print'Z_%d' % i
        os.chdir('Z_%d' % i)
        print 'In directory Z_%d' % i
        for filename1 in glob.glob('*.map'):
            os.remove(filename1)
            #removes the .map files
        print 'Removed maps'
        for filename2 in glob.glob('*.fld'):
            os.remove(filename2)
            #remove the .fld files
        print 'Removed fld'
        for filename3 in glob.glob('*.xyz'):
            os.remove(filename3)
            #removes the .xyz files
        print 'Removed xyz'
        os.chdir('../')
    else:
        print 'Z_%d directory does not exist' % i
        continue
print 'Done removing files generated to complete docking experiment'
#The .glg and .dlg files remain so that you can see what was used in generating grids
# and what was used during the docking experiment.
        