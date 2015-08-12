import MDAnalysis
import glob, os, string, sys
import subprocess
import shutil
from shutil import copy2
#sandbox_dock contains all the ligand PDBQTs
folder = '/PATH/TO/YOUR/sandbox_dock'

i = 0
j = 0
k = 0
m = 0
n = 0 
p = 0

#This loop takes the name of the file and makes it into a directory based on its basename.
# After the directory is created, the file corresponding to the directory name is moved into the corresponding directory.
for file_path in glob.glob(os.path.join(folder, '*.*')):
    new_dir = file_path.rsplit('.', 1)[0]
    os.mkdir(os.path.join(folder, new_dir))
    shutil.move(file_path, os.path.join(new_dir, os.path.basename(file_path)))
    i = i + 1
    print "Created directory number: ", i #Note that i does not correspond to the directory name, simply how many directories have been created.
print 'Done creating directories for docking, based on ligand name'

#set j to total number of directories created
for j in range(90):
    if os.path.exists('/PATH/TO/YOUR/sandbox_dock/NewDirectory%d' % j):
        os.chdir('/PATH/TO/YOUR/sandbox_dock/NewDirectory%d' % j)
        #this step prepares the ligand for docking. The -C flag preserves ligand charges; -l flag specifies the ligand to be prepared; -o flag specifies the output ligand filename. I recommend using something simple, like ligand.pdbqt.
        os.system('/PATH/TO/YOUR/MGLTools-1.5.6/bin/pythonsh /PATH/TO/YOUR/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l PDBQTFileStemName%d.pdbqt -C -o ligand.pdbqt' % j)
        #next step makes a symbolic link to your receptor PDB
        os.symlink('PATH/TO/PROTEIN/receptor.pdb', 'receptor.pdb')
        #This step prepares your receptor into the PDBQT format. The -A hydrogens add hydrogens to the protein, since this information generally lacks from crystallographic data most commonly used
        #to generate PDBs. The -r flag is to specify the receptor to be prepared and the -o flag specifies the name of the output file.
        os.system('/PATH/TO/YOUR/MGLTools-1.5.6/bin/pythonsh /PATH/TO/YOUR/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r receptor.pdb -A hydrogens -o receptor.pdbqt')
    else:
        print 'Directory NewDirectory%d does not exist!' % j
print 'Done linking ligands and receptor proteins'

# set k to total number of directories created
for k in range(90):
    if os.path.exists('/PATH/TO/YOUR/sandbox_dock/NewDirectory%d' % k):
        os.chdir('/PATH/TO/YOUR/sandbox_dock/NewDirectory%d' % k)
        #set m's range to the total number of docking experiments you wish to carry out
        #this step prepares the 'Config' file that AutoDock Vina uses in lieu of the grid maps, etc... to be generated. Adjust the parameters to fit your needs.
        #Check Vina documentation for additional options not specified below.
        for m in range(1,11):
            with open('receptor_%d_%d.txt' % (k, m), 'w') as fin:
                fin.write('receptor = receptor.pdbqt\nligand = ligand.pdbqt\ncenter_x = 1.445\ncenter_y = 7.314\ncenter_z = -7.982\nsize_x = 30.0\nsize_y = 22.5\nsize_z = 37.5\nout = receptor_%d_%d.pdbqt\nlog = receptor_%d_%d.log' % (k, m, k, m))
            #adjust these conditions so that the loop does indeed break, if not automatically.
            if m < 10 and k < 89:
                continue
            elif m == 10 and k == 89:
                break
    else:
        print 'NewDirectory%d does not exist!!' % k
#set n to total number of directories created
for n in range(90):
    if os.path.exists('/PATH/TO/YOUR/sandbox_dock/NewDirectory%d' % n):
        os.chdir('/PATH/TO/YOUR/sandbox_dock/NewDirectory%d' % n)
        #set p's range to the same range you did for m above, i.e. the total number of docking experiments to be carried out
        #--config flag is used to specify the config file with all the information necessary for docking to occur.
        #The following step will start your docking experiments. Vina automatically detects how many CPUs are available, and 
        #will use all that are available. If Vina cannot detect how many CPUs, it will use only 1 CPU for the computation.
        #Further development may want to carry out several docking experiements simultaneously.
        for p in range(1,11):
            os.system('/PATH/TO/YOUR/autodock_vina_1_1_2_linux_x86/bin/vina --config receptor_%d_%d.txt' % (n, p))
            #adjust these conditions so that the loop breaks, if it does not automatically do so already.
            if n < 10 and p < 89:
                print 'Still more dockings to go'
                continue
            elif n == 10 and p == 89:
                break
    else:
        print 'NewDirectory%d does not exist!!!' % n
print 'Docking experiments complete!'