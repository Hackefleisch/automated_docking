
from rdkit import Chem
from rdkit.Chem import AllChem
from ligand_diffusion import diffuse_ligand

from pyrosetta import *

pyrosetta.init()

def smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    return mol

def mol_to_res(mol):
    chem_manager = pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()

    tag = "fa_standard"
    lig_restype = pyrosetta.rosetta.core.chemical.MutableResidueType( 
        chem_manager.atom_type_set( tag ),
        chem_manager.element_set( "default" ),
        chem_manager.mm_atom_type_set( tag ),
        chem_manager.orbital_type_set( tag )
    )
    lig_restype.name( "UNKNOWN" )
    lig_restype.name3( "UNK" )
    lig_restype.name1( "X" )
    lig_restype.interchangeability_group( "UNK" )

    index_to_vd = {}

    conf = mol.GetConformer( 0 )

    for i in range( mol.GetNumAtoms() ):
        atom = mol.GetAtomWithIdx( i )
        element_name = atom.GetSymbol()
        charge = atom.GetFormalCharge()

        vd_atom = lig_restype.add_atom( "" )
        restype_atom = lig_restype.atom( vd_atom )
        restype_atom.element_type( lig_restype.element_set().element( element_name ) )
        restype_atom.formal_charge( charge )
        restype_atom.mm_name( "VIRT" )

        atom_pos = conf.GetAtomPosition( i )
        xyz = pyrosetta.rosetta.numeric.xyzVector_double_t( atom_pos.x, atom_pos.y, atom_pos.z )
        restype_atom.ideal_xyz( xyz )

        index_to_vd[ i ] = vd_atom


    for bond in mol.GetBonds():
        bond_name = bond.GetBondType()
        if bond_name == Chem.rdchem.BondType.SINGLE:
            bond_name = pyrosetta.rosetta.core.chemical.BondName.SingleBond
        elif bond_name == Chem.rdchem.BondType.DOUBLE:
            bond_name = pyrosetta.rosetta.core.chemical.BondName.DoubleBond
        elif bond_name == Chem.rdchem.BondType.TRIPLE:
            bond_name = pyrosetta.rosetta.core.chemical.BondName.TripleBond
        elif bond_name == Chem.rdchem.BondType.AROMATIC:
            bond_name = pyrosetta.rosetta.core.chemical.BondName.AromaticBond
        else:
            print( "ERROR: encountered unknown bond type", bond_name )
            bond_name = pyrosetta.rosetta.core.chemical.BondName.UnknownBond

        lig_restype.add_bond( 
            index_to_vd[ bond.GetBeginAtom().GetIdx() ],
            index_to_vd[ bond.GetEndAtom().GetIdx() ],
            bond_name
        )

    pyrosetta.rosetta.core.chemical.rename_atoms( lig_restype, True )
    pyrosetta.rosetta.core.chemical.rosetta_retype_fullatom( lig_restype, True )
    pyrosetta.rosetta.core.chemical.rosetta_recharge_fullatom( lig_restype )

    pyrosetta.rosetta.core.chemical.find_bonds_in_rings( lig_restype )

    nbr_vd = 0
    shortest_nbr_dist = 999999.99
    for vd in index_to_vd.values():
        if lig_restype.atom( vd ).element_type().get_chemical_symbol() == "H":
            continue
        tmp_dist = pyrosetta.rosetta.core.chemical.find_nbr_dist( lig_restype, vd )
        if tmp_dist < shortest_nbr_dist:
            shortest_nbr_dist = tmp_dist
            nbr_vd = vd

    lig_restype.nbr_radius( shortest_nbr_dist )
    lig_restype.nbr_atom( nbr_vd )
    lig_restype.assign_internal_coordinates()
    lig_restype.autodetermine_chi_bonds()

    lig_restype_non_mutable = pyrosetta.rosetta.core.chemical.ResidueType.make( lig_restype )
    return pyrosetta.rosetta.core.conformation.Residue( lig_restype_non_mutable, True )
    
def load_pose(path):

    return pyrosetta.pose_from_pdb( path )

def create_complex(pose, res):

    copy_pose = pyrosetta.rosetta.core.pose.Pose()
    copy_pose.detached_copy(pose)

    copy_pose.append_residue_by_jump( res, 1, "", "", True )
    # I will assume the last residue is the ligand
    copy_pose.pdb_info().chain( copy_pose.total_residue(), 'X' )
    copy_pose.update_pose_chains_from_pdb_chains()

    return copy_pose

def create_protocol(path, mover, scfx):

    xml_objects = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file(path)
    protocol = []
    for m in mover:
        protocol.append(xml_objects.get_mover(m))
    score_func = []
    for sf in scfx:
        score_func.append(xml_objects.get_score_function(sf))
    return protocol, score_func

class FullDocker:

    def __init__(self, smiles, pose, protocol, scfx) -> None:
        self.smiles = smiles
        self.pose = pose
        self.protocol = protocol
        self.scfx = scfx
        self.step = 0
        self.max_steps = 4
        self.mol = None
        self.distance = None
        self.next_step = "Diffusing ligand"
        self.max_num_dockings = 30
        self.finished_dockings = 0
        self.complex = None
        self.best_score = 9999.999
        self.best_complex = pyrosetta.rosetta.core.pose.Pose()
        self.work_pose = pyrosetta.rosetta.core.pose.Pose()

    # needs to be called several times until returning true
    def run(self) -> bool:
        if self.step == 0:
            self.mol, self.distance = diffuse_ligand(self.smiles)
            self.step += 1
            self.next_step = "Dock ligand step " + str(self.finished_dockings + 1) + "/" + str(self.max_num_dockings)
            return False
        elif self.step == 1:
            res = mol_to_res(self.mol)
            self.complex = create_complex(self.pose, res)
            self.best_complex.detached_copy(self.complex)
            self.step += 1
            return False
        elif self.step == 2:
            if self.finished_dockings < self.max_num_dockings:
                print("----Protocol round", self.finished_dockings + 1)
                self.work_pose.detached_copy(self.complex)
                for p in self.protocol:
                    print("----Apply", p.get_name())
                    p.apply(self.work_pose)
                idelta = pyrosetta.rosetta.protocols.ligand_docking.get_interface_deltas( 'X', self.work_pose, self.scfx[0] )
                score = idelta["interface_delta_X"] - idelta["if_X_coordinate_constraint"]
                if score < self.best_score:
                    self.best_score = score
                    self.best_complex.detached_copy(self.work_pose)
                print("----Current score:", score, "best score:", self.best_score)
                self.finished_dockings += 1
                self.next_step = "Dock ligand step " + str(self.finished_dockings + 1) + "/" + str(self.max_num_dockings)
                return False
            else:
                self.step += 1
                self.next_step = "Normalize"
                return False
        elif self.step == 3:
            n_atoms = Chem.rdMolDescriptors.CalcNumHeavyAtoms(self.mol)
            print("Number of heavy atoms:", n_atoms, "Normalization:", n_atoms ** 0.5, "old score:", self.best_score)
            self.best_score /= n_atoms ** 0.5
            print("new score", self.best_score)
            self.step += 1
            return True        


def full_docking(smiles, pose, protocol, scfx):
    mol, distance = diffuse_ligand(smiles)
    res = mol_to_res(mol)
    complex = create_complex(pose, res)
    best_score = 999999.99
    best_complex = pyrosetta.rosetta.core.pose.Pose()
    best_complex.detached_copy(complex)
    work_pose = pyrosetta.rosetta.core.pose.Pose()
    for pr in range(30):
        print("----Protocol round", pr)
        work_pose.detached_copy(complex)
        for p in protocol:
            print("----Apply", p.get_name())
            p.apply(work_pose)
        idelta = pyrosetta.rosetta.protocols.ligand_docking.get_interface_deltas( 'X', work_pose, scfx[0] )
        score = idelta["interface_delta_X"] - idelta["if_X_coordinate_constraint"]
        print("----Current score:", score, "best score:", best_score)
        if score < best_score:
            best_score = score
            best_complex.detached_copy(work_pose)
    n_atoms = Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol)
    best_score = best_score / (n_atoms**0.5)
    return best_complex, best_score, distance

def main():
    pose = load_pose( "input/8ef6_0001_R.pdb" )
    protocol, scfx = create_protocol("input/transform_repack.xml", ["high_res_docker", "final"], ["hard_rep"])

    smiles = "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
    complex, score, distancce = full_docking(smiles, pose, protocol, scfx)
    print(complex.scores)
    print(score)

    complex.dump_pdb("output/tmp.pdb")

if __name__ == "__main__":
    main()