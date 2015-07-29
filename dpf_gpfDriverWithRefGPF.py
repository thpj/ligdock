import os, string, sys
import subprocess

a = 0
b = 0
#The range for 'a' is determined by the number of directories. In our case this corresponded to the ligand ID number.
#The directories are setup such that they have the structure: dockingfiles/Z_%d so that it is clear which ligand we are docking, and for later analysis.
for a in range(30450):
    #change path to your directory with files for docking
    if os.path.exists('/Path/To/Your/Files/dockingfiles/Z_%d' % a):
        os.chdir('Z_%d' % a)
        #The range for 'b' is determined by how many types of receptors present. In this case, a range of 10 corresponds to receptor conformations labeled 0 to 9.
        for b in range(10):
            if os.path.exists('yourReceptorBasename%d.pdbqt' % b):
                #This step creates the GPF file. -l flag specifies the ligand; -r flag specifies the receptor; -i flag specifies the gpf file created before to use as reference; -o flag specifies the output file
                #This step assumes that you have created a reference gpf file using AutoDockTools. If you are targeting a specific region of the protein, this will be necessary; however, if blind docking you may just create a grid box that encompasses the entire protein instead.
                os.system('/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -l ligand.pdbqt -r yourReceptorBasename%d.pdbqt -i yourReceptorBasename%d.gpf -o yourReceptorBasename%d_%d.gpf' % (b, b, b, a))
                #This step creates the DPF file. -l flag specifies the ligand; -r flag specifies the recepto; -o specifies the output file
                os.system('/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -l ligand.pdbqt -r yourReceptorBasename%d.pdbqt -o yourReceptorBasename%d_%d.dpf' % (b, b, a))
                print 'Created yourReceptorBasename%d_%d gpf and dpf' % (b, a)
                if b == 9:
                    print 'Changing directories'
                    os.chdir('../')
                    print os.getcwd()
                else:
                    print 'Am here'
                    continue
            else:
                print 'yourReceptorBasename%d.pdbqt does not exist' % b
                continue
    else:
        print 'Directory Z_%d does not exist' % a
        continue
print 'Done creating yourReceptorBasename gpfs and dpfs'

print 'All grid and docking parameter files have been written!'