"""
Calculate the flexibility score from a PDB file of different origins.

Given a certain PDB file as input, the program identifies its source
(NMR, X-ray or AlphaFold prediction) and calculates an appropriate flexibility
score for each residue in described in the file (only considering backbone atoms).
"""
from cProfile import label
import copy
import numpy as np
from scipy.spatial.distance import squareform
from math import sqrt, pi, log10
import matplotlib.pyplot as plt


class Atom:
    """ Class storing the information given for each atom in the input PDB file. """

    def __init__(self, line):
        self.name = line[11:16].strip()
        self.resName = line[17:20]
        self.chainID = line[21:22]
        self.resSeq = int(line[22:26])
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.tempFactor = line[60:66].strip()
        if self.tempFactor:
            self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]


class PDB:
    """
    Object parsing the content of a PDB file from a protein obtained through NMR.

    It only stores information regarding lines starting with 'ATOM' as well as the
    number of model to which each set of atoms corresponds.
    """

    def __init__(self, file):
        self.file = file
        self.atoms = []
        self.type = None
        self.parse()

    def parse(self):
        """Store the characteristics of the NMR PDB file."""
        MODEL = None
        f = open(self.file, 'r')
        for line in f.readlines():
            if line.startswith('MODEL'):
                MODEL = int(line.split()[1])
            if line.startswith('EXPDTA'):
                self.type = str(line.split()[1])      
            if line.startswith('ATOM'):
                atom = Atom(line)
                atom.MODEL = MODEL
                self.atoms.append(atom)
        f.close()

    def get_atoms(self):
        """Return a list of all atoms in the PDB file"""
        return self.atoms


def fit_rms(ref_c, c):

    # Move geometric center of the structure to the origin of the coordinate
    # system
    ref_trans = np.average(ref_c, axis=0)
    ref_c = ref_c - ref_trans
    c_trans = np.average(c, axis=0)
    c = c - c_trans

    # Generate the covariance matrix
    C = np.dot(c.T, ref_c)

    # Obtain the Singular Value Decomposition
    (r1, s, r2) = np.linalg.svd(C)

    # Compute sign (remove mirroring)
    if np.linalg.det(C) < 0:
        r2[2, :] *= -1.0
    U = np.dot(r1, r2)
    return (c_trans, U, ref_trans)


class RMSDcalculator:
    """Object containing the methods necessary to calculate the RMSD value from a PDB NMR file."""

    def __init__(self, atoms1, atoms2, name=None):
        xyz1 = self.get_xyz(atoms1, name=name)
        xyz2 = self.get_xyz(atoms2, name=name)
        self.set_rmsd(xyz1, xyz2)

    def get_xyz(self, atoms, name=None):
        """Converts atom coordinates in the PDB file to an array object."""
        xyz = []
        for atom in atoms:
            if name:
                if atom.name != name:
                    continue
            xyz.append([atom.x, atom.y, atom.z])
        return np.array(xyz)

    def set_rmsd(self, c1, c2):
        """Calculates the RMSD value for two NMR structures."""
        self.rmsd = 0.0
        self.c_trans, self.U, self.ref_trans = fit_rms(c1, c2)
        new_c2 = np.dot(c2 - self.c_trans, self.U) + self.ref_trans
        self.rmsd = np.sqrt(np.average(np.sum((c1 - new_c2)**2, axis=1)))

    def get_aligned_coord(self, atoms, name=None):
        """Calculates the coordinates of a structure relative to a reference structure"""
        new_c2 = copy.deepcopy(atoms)
        for atom in new_c2:
            atom.x, atom.y, atom.z = np.dot(
                np.array([atom.x, atom.y, atom.z]) - self.c_trans, self.U) + self.ref_trans
        return new_c2


class NMR(object):
    """
    Object containing the parsing of NMR files.

    It calculates the average structure among the models and retrieves the model that has
    the highest similarity to it.
    """

    def __init__(self, atoms):
        # Read each MODEL
        self.MODELs = {}
        for atom in atoms:
            # Go through the different MODELs in the file
            if atom.MODEL not in self.MODELs:
                self.MODELs[atom.MODEL] = []
            # Choose only backbone atoms
            if atom.name not in ['C', 'CA', 'N', 'O']:
                continue
            self.MODELs[atom.MODEL].append(atom)

        # Set reference (using a centroid, qhich is one of the MODELs)
        ref_index = self.get_centroid()
        self.reference = self.MODELs[ref_index]
        self.set_RMSF(self.reference)

    def get_average_MODEL(self):
        """Calculates the average model among the coordinates of all the NMR models."""
        reference = self.MODELs[1]
        XYZ = []
        for n in self.MODELs:
            MODEL = self.MODELs[n]
            xyz = RMSDcalculator(reference, MODEL).get_aligned_coord(MODEL)
            xyz = np.array([[atom.x, atom.y, atom.z] for atom in xyz])
            XYZ.append(xyz)

        # Determine and return average structure
        m, n = xyz.shape
        sum_xyz = np.zeros((m, n))
        for xyz in XYZ:
            sum_xyz += xyz
        avg_xyz = sum_xyz/len(XYZ)
        avg_MODEL = copy.deepcopy(MODEL)
        for n in range(len(avg_MODEL)):
            avg_MODEL[n].x = avg_xyz[n][0]
            avg_MODEL[n].y = avg_xyz[n][1]
            avg_MODEL[n].z = avg_xyz[n][2]
        return avg_MODEL

    def get_centroid(self, avg_MODEL=None):
        """Return the MODEL number corresponding to the centroid."""
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

    def distx2(self, x1, x2):
        """Calculate the square of the distance between two coordinates."""
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
                # Differences between the state and the reference_state for each atom.
                diff[j] += self.distx2(np.array([self.MODELs[i][j].x,
                                                 self.MODELs[i][j].y, self.MODELs[i][j].z]), refer_xyz[j])
            diff[j] = np.sqrt(diff[j]/len(self.MODELs))

            if not (self.MODELs[i][j].chainID, self.MODELs[i][j].resSeq, self.MODELs[i][j].name) in tempFactor:
                tempFactor[(self.MODELs[i][j].chainID, self.MODELs[i]
                            [j].resSeq, self.MODELs[i][j].name)] = []
            tempFactor[(self.MODELs[i][j].chainID, self.MODELs[i][j].resSeq,
                        self.MODELs[i][j].name)].append(diff[j]*8*np.pi**2)

        RMSF = {}
        for n in range(natom):
            if not (reference[n].chainID, reference[n].resSeq, reference[n].resName) in RMSF:
                RMSF[(reference[n].chainID, reference[n].resSeq,
                      reference[n].resName)] = []
            RMSF[(reference[n].chainID, reference[n].resSeq,
                  reference[n].resName)].append(diff[n])
        RMSF = {k: np.mean(RMSF[k]) for k in RMSF}
        self.RMSF = RMSF

        self.tempFactor = tempFactor


def calculation_from_NMR(pdb_file, name_chain):
    """
    Calculate RMSF of each residue in a PDB file of a protein obtained by NMR.

    As a reference model for the RMSD computations, the centroid is used (the structure
    that is least different from the mean structure of all the models in the PDB file).
    """
    # Read PDB file
    pdb = PDB(pdb_file)
    # Get backbone atoms
    atoms = pdb.get_atoms()

    # Get the mean structure
    nmr = NMR(atoms)
    mean_structure = nmr.reference

    # Get the model number of the centroid structure
    centroid_model = nmr.get_centroid(mean_structure)

    # Write RMSF by residue to the output file
    RMSF_list = list(nmr.RMSF.items())

    # Sort the values by chain and residue number
    RMSF_list.sort()

    scores = []

    # Write the header of the output file
    scores.append(["ResName", "Chain", "ResID", "RMSF"])

    # Generate an array with flexibility scores
    flexibility = []
    number_of_residue = []
    for (chain, resID, resname), rmsf in RMSF_list:
        if chain == name_chain:
            flexibility.append(rmsf)
            number_of_residue.append(resID)

    # Normalise the flexibility scores of the array
    for (chain, resID, resname), rmsf in RMSF_list:
        if chain == name_chain:
            rmsf = (rmsf - min(flexibility)) / \
                    (max(flexibility)-min(flexibility))
            scores.append([resname, chain, resID, rmsf])
    return scores


def rmsf_Bfactor(atoms):
    """Returns a list of root-mean square fluctuations (RMSF) for each atom in the PDB."""

    # Extrac the B-factor of each atom in the PDB file
    atom_bfactors = [atom.tempFactor for atom in atoms]

    # Obtain the RMSF of each atom, with respect to the set
    def calc_rmsf(b): return sqrt(3 * b / (8 * pi**2))
    RMSFs = map(calc_rmsf, atom_bfactors)

    return RMSFs


def calculation_from_crystal(pdb_file, name_chain):
    """
    Calculate RMSF of each residue in a PDB file of a protein obtained by X-ray.

    The RMSF is calculated from the B-factor values (RMSF = sqrt(3*B/(8*pi**2))).
    """

    # Read PDB file
    pdb = PDB(pdb_file)

    # Extract only backbone atoms
    atoms = pdb.get_atoms()
    allowed_names = ['C', 'CA', 'N', 'O']
    backbone_atoms = [atom for atom in atoms if atom.name in allowed_names]

    # Calculate RMSF of all backbone atoms
    RMSF_list = rmsf_Bfactor(backbone_atoms)

    means, current_rmsfs = [], []
    counter = 0

    # Choose the adequate chain based on the matching PDB
    for atom in atoms:
        if (atom.chainID == name_chain) and (counter == 0):
            current_resid = atom.resSeq
            residue_list = [backbone_atoms[0]]
            counter = 1

    # Get a list of mean RMSFs for each residue
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        if atom.chainID == name_chain:
            if atom.resSeq != current_resid:
                means.append(sum(current_rmsfs) / len(current_rmsfs))
                residue_list.append(atom)
                current_rmsfs = []
                current_resid = atom.resSeq
            current_rmsfs.append(rmsf)
    means.append(sum(current_rmsfs) / len(current_rmsfs))
    RMSF_list = means
    backbone_atoms = residue_list
    keys = ['resName', 'chainID', 'resSeq']
    scores = []

    # Write the header of the output file
    scores.append(["ResName", "Chain", "ResID", "RMSF"])

    # Generate an array with flexibility scores
    flexibility = []
    number_of_residue = []
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        values = [atom[key] for key in keys]
        if values[1] == name_chain:
            flexibility.append(rmsf)
            number_of_residue.append(values[2])

    # Normalise the flexibility scores of the array
    new_flexibility = []
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        values = [atom[key] for key in keys]
        if values[1] == name_chain:
            rmsf = (rmsf - min(flexibility)) / \
                    (max(flexibility)-min(flexibility))
            new_flexibility.append(rmsf)
            values.append(rmsf)
            scores.append(values)
    return scores


def rmsf_pLLDT(atoms):
    """Returns a list of values assumed to RMSF for each atom in the PDB."""

    # Extrac the pLLDT value of each atom in the PDB file
    atom_pLLDT = [atom.tempFactor for atom in atoms]
    # Obtain the RMSF of each atom, with respect to the set
    def calc_rmsf(p): return log10(100-p)
    RMSFs = map(calc_rmsf, atom_pLLDT)
    return RMSFs


def calculation_from_alphafold(pdb_file):
    """
    Calculate RMSF of each residue in a PDB file predicted by AlphaFold.

    The RMSF is calculated from the pLDDT values (RMSF = log10(100-pLDDT)).
    """
    # Read PDB file
    pdb = PDB(pdb_file)

    # Extract only backbone atoms
    atoms = pdb.get_atoms()
    allowed_names = ['C', 'CA', 'N', 'O']
    backbone_atoms = [atom for atom in atoms if atom.name in allowed_names]

    # Calculate RMSF of all backbone atoms
    RMSF_list = rmsf_pLLDT(backbone_atoms)

    means, current_rmsfs = [], []
    current_resid = backbone_atoms[0].resSeq
    residue_list = [backbone_atoms[0]]
    # Get a list of mean RMSFs for each residue
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        if atom.resSeq != current_resid:
            means.append(sum(current_rmsfs) / len(current_rmsfs))
            residue_list.append(atom)
            current_rmsfs = []
            current_resid = atom.resSeq
        current_rmsfs.append(rmsf)
    means.append(sum(current_rmsfs) / len(current_rmsfs))
    RMSF_list = means
    backbone_atoms = residue_list
    keys = ['resName', 'chainID', 'resSeq']
    scores = []

    # Write the header of the output file
    scores.append(["ResName", "Chain", "ResID", "RMSF"])

    # Generate an array with flexibility scores
    flexibility = []
    number_of_residue = []
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        values = [atom[key] for key in keys]
        flexibility.append(rmsf)
        number_of_residue.append(values[2])

    # Normalise the flexibility scores of the array
    for atom, rmsf in zip(backbone_atoms, RMSF_list):
        values = [atom[key] for key in keys]
        rmsf = (rmsf - min(flexibility))/(max(flexibility)-min(flexibility))
        values.append(rmsf)
        scores.append(values)
    return scores


def general_calculation(pdb_file, name_chain = None):
    """Calls the functions necessary to obtain the flexibility score.

    Based on the mehod indicated in the PDB file to obtain the structure,
    it chooses the method that must be used to calculate the flexbility score.
    """
    # Read PDB file
    pdb = PDB(pdb_file)
    if pdb.type == "X-RAY" or pdb.type == "ELECTRON":
        matrix = calculation_from_crystal(pdb_file, name_chain)
    elif pdb.type == "SOLUTION":
        matrix = calculation_from_NMR(pdb_file, name_chain)
    else:
        matrix = calculation_from_alphafold(pdb_file)
    return matrix


def general_calculation_multiple(pdb_list, alphafold):
    """Calls the functions necessary to obtain the flexibility score of a list
    of PDB files.

    Based on the mehod indicated in the PDB file to obtain the structure,
    it chooses the method that must be used to calculate the flexbility score.
    """

    flexibility_array = np.full([len(pdb_list)+2, int(alphafold.end)+1], None)
    j = 0
    for element in pdb_list:
        print(f"Computing flexibility score on {element.identifier}...")
        matrix = general_calculation(element.path, element.chain)
        print("Done")
        flexibility_array[j][0]=element.identifier
        for i in range(1, len(matrix)-1):
            flexibility_array[j][matrix[i][2]] = matrix[i][3]
        j = j + 1
    flexibility_array[-2][0] = 'Mean'
    flexibility_array[-1][0] = 'std'
    for k in range(1,int(alphafold.end)+1):
        column = flexibility_array[:-2, k]
        column = list(filter(None, column))
        if column:
            flexibility_array[-2][k] = np.mean(column)
            flexibility_array[-1][k] = np.std(column)
        else:
            flexibility_array[-2][k] = None
            flexibility_array[-1][k] = None        
    return flexibility_array

def represent_data(matrix, out, pdb_matrix_alphafold = []):
    """
    Using a matrix with normalized flexibility scores it returns a graphical representation of them and a text file with the values. 
    Flexibility score (y-axis) vs. Reside number (x-asis). If it AlphaFold score and PDB score is provided
    it represents both in the same plot.
    """
    if len(pdb_matrix_alphafold) > 0:
        matrix = np.concatenate((np.asarray(matrix, dtype = "str").transpose(), pdb_matrix_alphafold)).transpose()
        mask = matrix[0:, -2] != None
        mask[0] = False
        plt.errorbar(list(map(float, matrix[mask, 2])), list(map(float, matrix[mask, -2])), yerr=list(map(float, matrix[mask, -1])), color='blue', fmt='.', ecolor = 'lightblue', elinewidth = 1, capsize=2, label='PDBs')
        plt.scatter(list(map(float, matrix[1:, 2])), list(map(float, matrix[1:, 3])), color='red', marker='.', label = 'AlphaFold')
        plt.legend(loc='lower left')
    else:
        matrix = np.array(matrix, dtype = "str")
        plt.scatter(list(map(float, matrix[1:, 2])), list(map(float, matrix[1:, 3])), color='red', marker='.')
    plt.title('RMSF Vs Residue number', fontsize=14)
    plt.xlabel('Residue number', fontsize=14)
    plt.ylabel('Normalised RMSF', fontsize=14)
    plt.grid(True)
    plt.savefig(out + ".png", dpi=600)
    plt.clf()    
    np.savetxt(out + '.out', matrix, fmt="%s" , delimiter="\t")
