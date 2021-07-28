from ChargeOptimization import * 
import numpy as np
import os,sys

dirtop="/home/sgr40/Classify/FUS4+CYS"
dirpar="/home/sgr40/Classify/FUS4+CYS/toppar/"
dirtar="/home/sgr40/Classify/FUS4+CYS/Wat-Int/"
psf = os.path.join(dirtop,"ligandrm.psf")
pdb = os.path.join(dirtop,"ligandrm_geom.pdb")
scalingfactor = 1.16
step=0.02;
# Get Lennard Jones Parameters (Charmm_36)
pars=["con.prm",
      "par_all36_carb.prm",
      "par_all36_cgenff.prm",
      "par_all36_lipid.prm",
      "par_all36m_prot.prm",
      "par_all36_na.prm",
      "par_interface.prm"]
# Target Data
target=["CON-ACC-C1.log",
        "CON-ACC-O-120a.log",
        "CON-ACC-O.log",
        "CON-DON-C1.log",
        "CON-DON-H2.log",
        "CON-DON-H4.log",
        "CON-DON-H6.log",
        "CON-DON-N.log",
        "CON-ACC-N.log",
        "CON-ACC-O-120b.log",
        "CON-ACC-S.log",
        "CON-DON-H1.log",
        "CON-DON-H3.log",
        "CON-DON-H5.log",
        "CON-DON-H7.log"]
# Here I will include the atoms charge constraints and weights to optimize 
atoms=[[0,-1,1,],[1,-1.0,1.0,],[1,-1.0,1.0,],[],[]]
TotalCharge = 0.0;
watsp="/home/sgr40/Classify/FUS4+CYS/Wat-Int/wat-sp.log"
consp="/home/sgr40/Classify/FUS4+CYS/Wat-Int/CON-sp-HF.log"
conmp2="/home/sgr40/Classify/FUS4+CYS/Wat-Int/CON-sp-MP2.log"

dirpars=[os.path.join(dirpar,parfil) for parfil in pars]
dirtarget=[os.path.join(dirtar,tarfil) for tarfil in target]

#Now we get the CON
resname="CON"
optimizer=OptimizeCharges(pdb=pdb,psf=psf,
                 params=dirpars,
                 wat_sp=watsp,
                 res_sp=consp,
                 res_mp2=conmp2,
                 target=dirtarget,
                 resname=resname,
                 atoms=[],
                 charge=0.0)

optimizer.ParseTopology();
# Here we set the data of the targets
optimizer.SetTargetData();
# Lennard contribution cte
optimizer.SetLennard();
optimizer.SetConstrainedCharge();
# Set Coulomb Cte
optimizer.SetCoulombConstant();
# Now make optimizatin
print(np.shape(optimizer.Jacobian()))
