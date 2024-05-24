
from rdkit import Chem
from rdkit.Chem import AllChem

from pyrosetta import *

pyrosetta.init()

def smiles_to_res(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
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

    pose.append_residue_by_jump( res, 1, "", "", True )
    # I will assume the last residue is the ligand
    pose.pdb_info().chain( pose.total_residue(), 'X' )
    pose.update_pose_chains_from_pdb_chains()

    return pose

def create_protocol(path):

    xml_objects = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file("input/transform_repack.xml")
    return xml_objects.get_mover('ParsedProtocol')