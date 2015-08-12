import MDAnalysis
import os, string, sys
import subprocess
import shutil
import mdtraj as md

from MDAnalysis import *
from shutil import copy2

#This first part uses mdtraj to save a new DCD file, starting at frame 1 and skipping by 10 therafter. 
#This allows for the frames to be sequentially numbered for later use with MDAnalysis.
#This step may be skipped depending on your needs and resources available.
t = md.load('/PATH/TO/YOUR/fileName.dcd', top='/PATH/TO/YOUR/topology.pdb')
t[::10].save_dcd('/PATH/TO/YOUR/newFileName.dcd')
# save new trajectory with 10% of frames

u = MDAnalysis.Universe("/PATH/TO/YOUR/topology.pdb", "/PATH/TO/YOUR/newFileName.dcd", permissive=False)
#u loads your simulation into the MDAnalysis universe
print "Done Loading"
#iterates over trajectory; at this step, you can also specify at this point which frame you would like to start and stop at, and if you
#need to skip by a certain number of frames. Example: [start:stop:skip]
for ts in u.trajectory:
    ts = ts.frame
    # get the timestep of current frame
    os.system("mkdir Z_" + str(ts))
    # make directory with current timestep as string
    pdb = MDAnalysis.coordinates.PDB.PDBWriter("/PATH/TO/YOUR/Z_%d/PDBFileName_%d.pdb" % (ts, ts), universe=u, start=0, multiframe=True)
    #grabs the coordinates of the current timestep
    pdb.write(u)
    # setup PDB writer and write current frame as a PDB
    os.symlink("/PATH/TO/YOUR/Z_%d/PDBFileName%d.pdb" % (ts, ts), "/PATH/TO/YOUR/Z_%d/receptor.pdb" % (ts))
    # link the created PDB as receptor.pdb
    shutil.copy2("/PATH/TO/YOUR/ligandFileName.pdbqt", "/PATH/TO/YOUR/Z_%d/ligandFileName.pdbqt" % (ts))
    # copy the ligand to the newly created directory
    # if you have many ligands, it may be best to create links rather than a copy
    os.system("/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r /PATH/TO/YOUR/Z_%d/receptor.pdb -A hydrogens -o /PATH/TO/YOUR/Z_%d/receptor.pdbqt" % (ts, ts))
    # prepare the receptor PDBQT; -r specifies the receptor PDB, -A hydrogens adds hydrogens not present in the crystallographic informations, and -o specifies the output filename.
    os.system("/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -l /PATH/TO/YOUR/Z_%d/ligandFileName.pdbqt -r /PATH/TO/YOUR/Z_%d/receptor.pdbqt -p npts='160,160,160' -o /PATH/TO/YOUR/Z_%d/gridParameterFileName.gpf" % (ts, ts, ts))
    # prepare the grid parameter files of the receptor, given the ligand probe atom types. -l flag specifies the ligand, -r flag specifies the receptor, -p flage allows you to change parameters, such as the number of points in the x,y,and z directions, and the -o specifies the output filename.
    os.system("/PATH/TO/MGLTools-1.5.6/bin/pythonsh /PATH/TO/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -l /PATH/TO/YOUR/Z_%d/ligandFileName.pdbqt -r /PATH/TO/YOUR/Z_%d/receptor.pdbqt -o /PATH/TO/YOUR/Z_%d/dockingParameterFileName.dpf" % (ts, ts, ts))
    # prepare docking parameter files. -l flag specifies the ligand, -r flag specifies the receptor, and -o flag specifies the output file.
pdb.close()
print "All done"  
