import numpy as np
import re
import pickle 
import sys,os

class Atom:
    index=None;
    name=None;
    resname=None;
    x=None;
    y=None;
    z=None;
    el=None;
    mass=None;
    charge=None;
    atomtype=None;
    eps=None;
    r2=None;
    
class Target:
    refsep=None;
    distance_type=[];
    dihedrals=[];
    distances=[];
    dipole_exp=[];
    energy=[];
    
class OptimizeCharges:
    def __init__(self,pdb,psf=None,
                 params=[],wat_sp=None,res_sp=None,
                 res_mp2=None,target=[],resname=None,
                 atoms=[],charge=0.0,optstep=0.05):
        if psf==None:
            print("Reading ",pdb)
            # if psf is None then read the pdb as input fil
            # read the input file
        else:
            self.pdbf=pdb;
            self.psff=psf;
            self.paramsf=params;
            self.targetf=target;
            self.wat_spf=wat_sp;
            self.res_spf=res_sp;
            self.res_mp2f=res_mp2;
            self.resname=resname;
            self.atoms=atoms;
            self.tot_charge=charge;
            self.optstep=optstep;
            # Check if input is fine
            for i in [[self.wat_spf,"Water Single Point"],
                      [self.res_spf,"Residue Single Point"],
                      [self.res_mp2f,"Residue MP2"],[self.resname,"Resname"]]:
                if i[0]==None: 
                    print("The file for "+i[1]+" has not been specified.")
                    print("The optimization will not run without these files.")
            if len(self.paramsf)==0:
                print("The files for the target data have not been specified")
                print("The optimization will not run without these files.")
            if len(self.targetf)==0:
                print("All the charges will be used for the optimization")
        
        # End Check and Initialize variables that will be used during the Optimization 
        self.coordtop = None;
        self.refenergy= None;
        self.con_charge = None;
        self.coul_cte = None;
        self.atoms_constrained=[];
        self.qmint=[];
        self.watops=[];
        self.lennard=[];
        
        
        
    def RunOpt(self):
        pass
    
    def SetConstrainedCharge(self):
        if len(self.atoms)==0:
            print("All the atoms will be optimized in an constrained way between -1 and 1 and weight 1.")
            # Change the indexing
            self.con_charge=self.tot_charge;
            self.atoms=[[c,-1.0,1.0,atom.charge,1.0] for c,atom in enumerate(self.coordtop)]
            
        else: 
            # Here we will give the proper indication of the charges to optimize 
            indexes = [atom.index for atom in self.coordtop]
            constrained_atoms= [c for c,atom in enumerate(self.atoms) if atom not in indexes]
            to_optimize = [c for c,atom in enumerate(self.atoms) if atom in indexes]
            # Change the indexing 
            for c,i in enumerate(self.atoms):
                i[0]=to_optimize[c];
            # Get the constrained data
            constrained_charges = [self.coordtop[index].charge  for index in constrained_atoms]
            self.con_charge=self.tot_charge-np.sum(constrained_charges);
        return 1;
    
    def ParseTopology(self):
        # Read PDB
        st=self.ReadCoord();
        # Parse PSF
        if st==1:
            st=self.ParseCharges();
        else:
            print("Problems reading the Topology.")
            return 0;
        # Parse PAR
        if st==1:
            st=self.ParseLennard();
        else:
            print("Problems reading the Parameters.")
            return 0;
        
        return 1;

    def SetTargetData(self):
        # Here we set the reference energy. 
        self.SetRefEnergy();
        # Read the Target Data
        target_data = [self.Read_OPT_Gaussian(fil) for fil in self.targetf];
        # Water topologies and positions
        self.watops=[i[1] for i in target_data];
        # QM interaction energies
        self.qmint = np.array([(i[0]-self.refenergy)*627.50  for i in target_data])
        return 1;
    
    def SetRefEnergy(self):
        watenergy=self.Read_SP_Gaussian(self.wat_spf)
        conenergy=self.Read_SP_Gaussian(self.res_spf)
        self.refenergy=watenergy+conenergy;      
        return 1
    
    def SetLennard(self):
        # Here we compute the lennard jones cotribution for each water molecule.
        self.lennard = np.array([self.Lennard(self.coordtop,wat) for wat in self.watops])
        return 1
    
    def SetCoulombConstant(self):
        self.coul_cte=self.CoulombConstant();
        return 1
    
            
    def ReadCoord(self):
        atoms=[]
        fil=open(self.pdbf,"r")
        lines=fil.readlines()
        fil.close()
        for line in lines:
            if line[0:4]=="ATOM" and line[17:20].strip()==self.resname:
                atom=Atom()
                # need the element to detect the charge
                atom.index=int(line[6:11]);
                atom.name=line[12:16].strip();
                atom.resname=line[17:20].strip();
                atom.x=float(line[30:38]);
                atom.y=float(line[38:46]);
                atom.z=float(line[47:54]);
                atoms.append(atom)
            else:
                #Not do anything
                continue;
        #Update coordtop
        self.coordtop=atoms;
        return 1;

    
    def ParseCharges(self):
        # We do not need the topology as the structure is fixed
        mode=0;
        fil=open(self.psff,"r")
        lines=fil.readlines()
        fil.close()
        for line in lines:
            if mode==0:
                if line.find("!NATOM")!=-1:
                    start_atom=0;
                    nlines=int(line[0:11]);
                    count=0;
                    mode=1;
                    continue;
                pass
            elif mode==1:
                # If count get bigger than nlines stop process
                if count>=nlines:
                    mode=0;
                    continue;
                # If the line index coincides with the atom index the assign charge and mass
                # atomtypes
                print(line[0:11])
                if int(line[0:11])==self.coordtop[start_atom].index:
                    self.coordtop[start_atom].charge = float(line[56:71])
                    self.coordtop[start_atom].mass = float(line[72:79])
                    self.coordtop[start_atom].atomtype=line[47:53].strip()
                    # increase the atom index
                    start_atom=start_atom+1;
                count=count+1;
                pass
            else:
                continue
        return 1;
    
    def ParseLennard(self):
        # Get the unique 
        alltypes=np.array([at.atomtype for at in self.coordtop])
        uniquetypes = np.unique(alltypes)
        # READ ALL THE LENNARD JONES PARAMETERS TO A DICTIONARY
        dic_lennard={}
        for parfil in self.paramsf:
            mode=0
            fil=open(parfil,"r");
            lines=fil.readlines();
            fil.close()
            # Find in the file the keywords ATOMS BONDS ANGLES DIHEDRALS IMPROPERS NONBONDED NBFIX END
            keywords=["ATOMS","BONDS","ANGLES","DIHEDRALS","IMPROPERS","NONBONDED","NBFIX","END"]
            dic_how={}
            count=0
            for line in lines:
                for key in keywords:
                    if line[0:len(key)]==key:
                        dic_how[key]=count;
                count=count+1
            #Just read the NONBONDED
            getkeys = np.array(list(dic_how.keys()));
            index=[c1 for c1,t in enumerate(getkeys) if t=="NONBONDED"]
            if len(index)==1:
                index=index[0]
            else:
                continue
            start=dic_how["NONBONDED"]+2
            finish=dic_how[getkeys[index+1]]
            for line in lines[start:finish]:
                if len(line.strip())==0: continue
                if line.strip()[0]=="!": continue
                else:
                    table=re.split(" +",line)
                    atomtype=table[0]
                    eps=float(table[2]);
                    r2=float(table[3]);
                    if atomtype not in list(dic_lennard.keys()):
                        dic_lennard[atomtype]=[eps,r2]
                    else:
                        print("Twice lennard potential for atomtype",atomtype)
                        dic_lennard[atomtype]=[eps,r2]
            print("Finish Reading",parfil)
        # Assign After Reading All Parameter Files
        for at_un in uniquetypes:
            if at_un in list(dic_lennard.keys()):
                len_par=dic_lennard[at_un]
                logical=alltypes==at_un
                for c,log in enumerate(logical):
                    if log: 
                        self.coordtop[c].eps=len_par[0];
                        self.coordtop[c].r2=len_par[1];
            else:
                print("Missing Lennard parameter for atom: ",at_un)
                print("Add the missing in other Parameter Files")
                print("Returning empty list")
                # Something is missing
                return 0;
        # Everything went fine
        return 1;

    def Read_SP_Gaussian(self,logfile):
        fil = open(logfile,"r");
        lines = fil.readlines();
        for line in lines:
            if line.find("SCF Done")!=-1:
                energy=float(line[21:37].strip())
                break;
        # get the energy for the file
        return energy;

    def Read_OPT_Gaussian(self,logfile):
        # Get the optimized geometries
        fil=open(logfile,"r");
        energies=[];
        opt_orientation=[]
        lines=fil.readlines();
        fil.close()
        mode=0;
        # Read first all the energies
        for c,line in enumerate(lines):
            # Get also the starting structure and convert it with the Z matrix to some coordinates 
            if line.find("SCF Done")!=-1:
                energy=float(line[21:37].strip());
                energies.append(energy);
                continue;
            # Get the possitions of the minimized waters
            if line.find("Z-Matrix orientation")!=-1:
                atoms=[]
                for l in  lines[c+5:]:
                    if l.strip()[0]=="-":
                        break;
                    else:
                        table=re.split(" +",l.strip())
                        atoms.append(table)
                # Take the Last Four Positions and Append
                water_pos=atoms[-4:]
                opt_orientation.append(water_pos)
        # Get the last of the Optimization Steps and the Last energy
        opt=opt_orientation[-1]
        water=[]

        for atom in opt:
            newatom=Atom()
            newatom.resname="TIP3P"
            #print(atom[3],atom[4],atom[5])
            newatom.x= float(atom[3]);
            newatom.y= float(atom[4]);
            newatom.z= float(atom[5]);
            if int(atom[1])==-1:
                continue;
            elif int(atom[1])==8:
                newatom.atomtype="OT"
                newatom.name="O";
                newatom.charge=-0.834;
                newatom.eps=-0.1521;
                newatom.r2=1.7682;
                water.append(newatom);        
            elif int(atom[1])==1:
                newatom.atomtype="HT"
                newatom.name="H";
                newatom.charge=0.417;
                newatom.eps=-0.046;
                newatom.r2=0.2245;
                water.append(newatom);
            else:
                print("An atom of water is not H or O in ",logfile)
                print("Problem reading the the logfile ",logfile)
                sys.exit();
        energy=energies[-1]
        #print(energies[-1])
        print("The logfile readed correctly ",logfile)
        return (energy,water)

    def Lennard(self,comp_coortop,water_coortop):
        # This gives you the Lennard-Jones contribution to
        # The interaction energy with the water 
        # It is cte in time (execute just once)
        lenn=0;
        for i in comp_coortop:
            for j in water_coortop:
                ra = np.array([i.x,i.y,i.z]);
                rw = np.array([j.x,j.y,j.z]);
                d = np.linalg.norm(ra-rw)
                # Bertholot M. Rule
                eps  = np.sqrt(i.eps*j.eps);
                # Lorentz M. Rule (convert r/2 to sigma 1.78)
                sigma = (1.7817974/2.0)*(i.r2+j.r2);
                dsix = (sigma/d)**6
                lenn=lenn + eps * (dsix**2-dsix)

        lenn=4.0*lenn;
        return lenn;
    
    def CoulombConstant(self):
        if len(self.atoms_constrained)==0:
            return np.zeros(len(self.watops));
        else:
            # This will give a vector (a value for each target data)
            coul=0
            conv_factor=(627.50/14.3996)
            whole=[];
            for watop in self.watops:
                sum_count=0
                for i in self.atoms_constrained:
                    for j in watop:
                        atom=self.coordtop[i]
                        # Coulomb
                        ra = np.array([atom.x,atom.y,atom.z]);
                        rw = np.array([j.x,j.y,j.z]);
                        d = np.linalg.norm(ra-rw)
                        coul=coul+(atom.charge*j.charge/d)
                        sum_count=sum_count+coul       
                whole.append(sum_count)
            whole=conv_factor*np.array(whole)
            return whole     
        
    def CoulombMatrix(self):
        # This will give as a matrix (a value for each atom charged
        coul=0;
        conv_factor=0.332063714;
        whole=[]
        for watop in self.watops:
            mat=[]
            for i in self.atoms:
                atom=self.coordtop[i[0]]
                sum_count=0;
                for j in watop:
                    # Coulomb
                    ra = np.array([atom.x,atom.y,atom.z]);
                    rw = np.array([j.x,j.y,j.z]);
                    d = np.linalg.norm(ra-rw)
                    coul=coul+(i[3]*j.charge/d)
                    sum_count=sum_count+coul
                mat.append(sum_count)
            whole.append(np.array(mat))
        whole = conv_factor*np.array(whole);
        return whole;
    
    def MMEnergy(self,cmat):
        return np.sum(self.Cou,axis=1)+self.coul_cte+self.lennard;
        
    def IntEnergyCont(self):
        energy_cont=np.abs(self.qmint-self.MMEnergy());
        return energy_cont;
    
    def CoulombSingleAtom(self,atom,charge,watop):
        coul=0;
        conv_factor=0.332063714;
        for j in watop:
            ra = np.array([atom.x,atom.y,atom.z]);
            rw = np.array([j.x,j.y,j.z]);
            d = np.linalg.norm(ra-rw)
            coul=coul+(charge*j.charge/d)
        coul=conv_factor*coul;
        return coul;

    
    def Jacobian(self):
        coulmat = self.CoulombMatrix();
        Fx1=[]
        for watop in self.watops:
            inter=[];
            for atom in self.atoms:
                coords = self.coordtop[atom[0]]
                charge = atom[3]+self.optstep;
                inter.append(self.CoulombSingleAtom(coords,charge,watop))
            Fx1.append(inter);
        Fx1=np.array(Fx1);
        print(np.shape(coulmat))
        print(np.shape(Fx1))
        return Fx1; 
            # Get the row for the jacobian matrix. 
            

        # get the CoulombMatrix with charge+step
        

        #newcharge = charge0 - h*JACOBIAN  
        pass
