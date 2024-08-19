from ash import *

#Original raw PDB-file (no hydrogens, nosolvent) and xml
pdbfile = r"/scratch/09069/dhp563/ash/FeS_Az_nocluster/Az_CIACGAC_nocluster.pdb"

# Setting up system via Modeller
openmmobject, ashfragment = OpenMM_Modeller(pdbfile = pdbfile, 
                                            forcefield = "CHARMM36", 
                                            watermodel = "tip3p",
                                            pH = 7.0,
                                            solvent_padding = 10.0,
                                            ionicstrength = 0.1) 


#MM minimization for 1000 steps
OpenMM_Opt(fragment=ashfragment, theory=openmmobject, maxiter=1000, tolerance=1)

#Classical MD simulation for 1000 ps
OpenMM_MD(fragment=ashfragment, 
theory=openmmobject, 
timestep=0.001, 
simulation_time=1000, 
traj_frequency=1000, 
temperature=300,
integrator='LangevinMiddleIntegrator', 
coupling_frequency=1, 
trajectory_file_option='DCD')

