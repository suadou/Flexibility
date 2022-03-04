"""
Calculate the flexibility score from a PDB file of different origins
"""
import copy
import numpy as np
from scipy.spatial.distance import squareform
from math import sqrt, pi, log10
import matplotlib.pyplot as plt

class Atom:
    def __init__(self, line):
        self.name = line[11:16].strip()
        self.resName = line[17:20]
        self.chainID = line[21:22]
        self.resSeq = int(line[22:26])
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.tempFactor = line[60:66].strip()
        if self.tempFactor: self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]

class PDB:
    def __init__(self, file):
        self.file = file
        self.atoms = []
        self.parse()

    def parse(self):
        MODEL = None
        f = open(self.file, 'r')
        for line in f.readlines():
            if line.startswith('MODEL'): MODEL = int(line.split()[1])
            if line.startswith('ATOM'):
                atom = Atom(line)
                atom.MODEL = MODEL
                self.atoms.append(atom)
        f.close()

    def get_atoms(self, to_dict=True):
        """Return a list of all atoms.

        If to_dict is True, each atom is represented as a dictionary.
        Otherwise, a list of Atom objects is returned."""
        if to_dict: return [x.__dict__ for x in self.atoms]
        else: return self.atoms

    def get_model(self, model_num, to_dict=True):
        """Return all atoms where MODEL == model_num"""
        model_atoms = [x for x in self.atoms if x.MODEL == model_num]
        if to_dict:
            return [atom.__dict__ for atom in model_atoms]
        else:
            return model_atoms

class RMSDcalculator:
    def __init__(self, atoms1, atoms2, name=None):
        xyz1 = self.get_xyz(atoms1, name=name)
        xyz2 = self.get_xyz(atoms2, name=name)
        self.set_rmsd(xyz1, xyz2)

    def get_xyz(self, atoms, name=None):
        xyz = []
        for atom in atoms:
            if name:
                if atom.name != name: continue
            xyz.append([atom.x, atom.y, atom.z])
        return np.array(xyz)

    def set_rmsd(self, c1, c2):
        self.rmsd = 0.0
        self.c_trans, self.U, self.ref_trans = fit_rms(c1, c2)
        new_c2 = np.dot(c2 - self.c_trans, self.U) + self.ref_trans
        self.rmsd = np.sqrt( np.average( np.sum( ( c1 - new_c2 )**2, axis=1 ) ) )

    def get_aligned_coord(self, atoms, name=None):
        new_c2 = copy.deepcopy(atoms)
        for atom in new_c2:
            atom.x, atom.y, atom.z = np.dot(np.array([atom.x, atom.y, atom.z]) - self.c_trans, self.U) + self.ref_trans
        return new_c2

class NMR(object):
    def __init__(self, atoms):
        #read each MODEL
        self.MODELs = {}
        for atom in atoms:
            if atom.MODEL not in self.MODELs:
                self.MODELs[atom.MODEL] = []
            if atom.name not in ['C', 'CA', 'N', 'O']:
                continue
            self.MODELs[atom.MODEL].append(atom)

        # set reference(using centroid (one of the MODELs))

        ref_index = self.get_centroid()
        self.reference = self.MODELs[ref_index]


        self.set_RMSF(self.reference)

    def get_average_MODEL(self):
        reference = self.MODELs[1]
        XYZ = []
        for n in self.MODELs:
            MODEL = self.MODELs[n]
            xyz = RMSDcalculator(reference, MODEL).get_aligned_coord(MODEL)
            xyz = np.array([[atom.x, atom.y, atom.z] for atom in xyz])
            XYZ.append(xyz)

        # set average structure
        m,n = xyz.shape
        sum_xyz = np.zeros((m,n))
        for xyz in XYZ:
            sum_xyz += xyz
        avg_xyz = sum_xyz/len(XYZ)
        avg_MODEL = copy.deepcopy(MODEL)
        for n in range(len(avg_MODEL)):
            avg_MODEL[n].x = avg_xyz[n][0]
            avg_MODEL[n].y = avg_xyz[n][1]
            avg_MODEL[n].z = avg_xyz[n][2]
        return avg_MODEL

    def get_centroid(self, avg_MODEL = None):
        """Return the model number corresponding to the centroid."""
        centroid = None
        min_RMSD = None
        
        if not avg_MODEL:
            avg_MODEL = self.get_average_MODEL()
        for n in self.MODELs:
            rmsd = RMSDcalculator(avg_MODEL, self.MODELs[n]).rmsd
            if min_RMSD is None or rmsd < min_RMSD:
                n_centroid = n
                min_RMSD = rmsd
        return n_centroid

    def distx2(self,x1,x2):
        """Calculate the square of the distance between two coordinates.
        
        Returns a float."""
        distx2 = (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2
        return distx2
    
    def set_RMSF(self, reference):
        """Sets the self.RMSF attribute."""
        refer_xyz = np.array([[atom.x, atom.y, atom.z] for atom in reference])
        natom = len(refer_xyz)

        diff = []
        tempFactor = {}
        for j in range(natom):
            diff.append(0.)
            for i in self.MODELs:
                # differences between the states and reference_state for each atom.
                diff[j] += self.distx2(np.array([self.MODELs[i][j].x, self.MODELs[i][j].y, self.MODELs[i][j].z]), refer_xyz[j]) # d^2
            diff[j] = np.sqrt(diff[j]/len(self.MODELs))

            if not (self.MODELs[i][j].chainID, self.MODELs[i][j].resSeq, self.MODELs[i][j].name) in tempFactor:
                tempFactor[(self.MODELs[i][j].chainID, self.MODELs[i][j].resSeq, self.MODELs[i][j].name)] = []
            tempFactor[(self.MODELs[i][j].chainID, self.MODELs[i][j].resSeq, self.MODELs[i][j].name)].append(diff[j]*8*np.pi**2)

        RMSF = {}
        for n in range(natom):
            if not (reference[n].chainID, reference[n].resSeq, reference[n].resName) in RMSF:
                RMSF[(reference[n].chainID, reference[n].resSeq, reference[n].resName)] = []
            RMSF[(reference[n].chainID, reference[n].resSeq, reference[n].resName)].append(diff[n])
        RMSF = {k:np.mean(RMSF[k]) for k in RMSF}
        self.RMSF = RMSF
        
        self.tempFactor = tempFactor

def fit_rms(ref_c,c):
    # move geometric center to the origin
    ref_trans = np.average(ref_c, axis=0)
    ref_c = ref_c - ref_trans
    c_trans = np.average(c, axis=0)
    c = c - c_trans

    # covariance matrix
    C = np.dot(c.T, ref_c)

    # Singular Value Decomposition
    (r1, s, r2) = np.linalg.svd(C)

    # compute sign (remove mirroring)
    if np.linalg.det(C) < 0:
        r2[2,:] *= -1.0
    U = np.dot(r1, r2)
    return (c_trans, U, ref_trans)

def calculation_from_NMR(pdb_file):
    """
    Calculate RMSF of each residue in a PDB file of a protein obtained by NMR. As a reference model for the RMSD computations, the centroid is used (the structure that is least different from the mean structure of all the models in the PDB file)
    """
    # read PDB file
    pdb = PDB(pdb_file)
    # get backbone atoms
    atoms = pdb.get_atoms(to_dict=False)

    # get the mean structure
    nmr = NMR(atoms)
    mean_structure = nmr.reference

    # get the model number of the centroid structure
    centroid_model = nmr.get_centroid(mean_structure)
    ##print("# model", centroid_model, "is the representative structure")

    # display RMSF by residue
    RMSF_list = list(nmr.RMSF.items())
    RMSF_list.sort()  # sort by chain and residue number
    print("ResName\tChain\tResID\tRMSF")
    template = '{:>3s} {:>2s} {:>3d}  {:<f}'
    #number_of_residue = []
    #flexibility = []
    #for (chain, resID, resname), rmsf in RMSF_list:
    #    print(template.format(resname, chain, resID, rmsf))
    #    number_of_residue.append(resID)
    #    flexibility.append(log10(rmsf))
    #plt.plot(number_of_residue, flexibility, color='red', marker='o')
    #plt.title('RMSF Vs Residue number', fontsize=14)
    #plt.xlabel('Residue number', fontsize=14)
    #plt.ylabel('log10(RMSF)', fontsize=14)
    #plt.grid(True)
    #plt.show()

        

def rmsf_Bfactor(atoms):
    """Returns a list of root-mean square fluctuations.

    By default, the RMSF of each atom in atoms is returned. If by_residue
    is True, then the average RMSF for each residue is returned instead."""
    # B-factor of each atom
    atom_bfactors = [atom.tempFactor for atom in atoms]
    
    # RMSF of each atom, with respect to the set
    calc_rmsf = lambda b: sqrt(3 * b / (8 * pi**2))
    RMSFs = map(calc_rmsf, atom_bfactors)

    return RMSFs
    
def calculation_from_crystal(pdb_file):
    """
    Calculate the RMSF value of each residue in a PDB file with data from the B-factor (RMSF = sqrt(3*B/(8*pi**2)))
    """
    # read PDB file
    pdb = PDB(pdb_file)
    # get backbone atoms
    atoms = pdb.get_atoms(to_dict=False)
    allowed_names = ['C', 'CA', 'N', 'O']
    backbone_atoms = [atom for atom in atoms if atom.name in allowed_names]

    # get RMSF of all backbone atoms
    RMSF_list = rmsf_Bfactor(backbone_atoms)

    means, current_rmsfs = [], []
    current_resid = backbone_atoms[0].resSeq
    residue_list = [backbone_atoms[0]]
    # get a list of mean RMSFs for each residue
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        if atom.resSeq != current_resid:
            means.append(sum(current_rmsfs) / len(current_rmsfs))
            residue_list.append(atom)
            current_rmsfs = []
            current_resid = atom.resSeq
        current_rmsfs.append(rmsf)
    RMSF_list = means
    backbone_atoms = residue_list
    keys = ['resName', 'chainID', 'resSeq']

    print("ResName\tChain\tResID\tRMSF")

    # template string for justified output columns
    template = '{:>3s} {:>2s} {:>3d}  {:<f}'

    # output RMSF data for all backbone atoms
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        values = [atom[key] for key in keys]
        values.append(rmsf)
        print(template.format(*values))

def rmsf_pLLDT(atoms):
    """Returns a list of root-mean square fluctuations.

    By default, the RMSF of each atom in atoms is returned. If by_residue
    is True, then the average RMSF for each residue is returned instead."""
    # B-factor of each atom
    atom_pLLDT = [atom.tempFactor for atom in atoms]

    # RMSF of each atom, with respect to the set
    calc_rmsf = lambda p: log10(p)
    RMSFs = map(calc_rmsf, atom_pLLDT)

    return RMSFs

def calculation_from_alphafold(pdb_file):
    """
    Calculate the log10(pLDDT) of each residue of a PDB file obtained from AlphaFold
    """
    # read PDB file
    pdb = PDB(pdb_file)
    # get backbone atoms
    atoms = pdb.get_atoms(to_dict=False)
    allowed_names = ['C', 'CA', 'N', 'O']
    backbone_atoms = [atom for atom in atoms if atom.name in allowed_names]

    # get RMSF of all backbone atoms
    RMSF_list = rmsf_pLLDT(backbone_atoms)

    means, current_rmsfs = [], []
    current_resid = backbone_atoms[0].resSeq
    residue_list = [backbone_atoms[0]]
    # get a list of mean RMSFs for each residue
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        if atom.resSeq != current_resid:
            means.append(sum(current_rmsfs) / len(current_rmsfs))
            residue_list.append(atom)
            current_rmsfs = []
            current_resid = atom.resSeq
        current_rmsfs.append(rmsf)
    RMSF_list = means
    backbone_atoms = residue_list
    keys = ['resName', 'chainID', 'resSeq']

    print("ResName\tChain\tResID\tRMSF")

    # template string for justified output columns
    template = '{:>3s} {:>2s} {:>3d}  {:<f}'

    # output RMSF data for all backbone atoms
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        values = [atom[key] for key in keys]
        values.append(rmsf)
        print(template.format(*values))
