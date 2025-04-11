from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from pyrosetta import *
from openbabel import pybel
from matplotlib import pyplot as plt
from rdkit.Chem import Descriptors
from rdkit.Chem import QED

pyrosetta.init()

def rdkit_to_mutable_res(mol):
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

    index_to_name = {}

    for idx in range( mol.GetNumAtoms() ):
        vd = index_to_vd[idx]
        atm_name = lig_restype.atom_name( vd )
        index_to_name[ idx ] = atm_name

    return lig_restype, index_to_vd, index_to_name

def add_confs_to_res(mol, mutable_res, index_to_vd):

    rotamers_spec = rosetta.core.chemical.rotamers.StoredRotamerLibrarySpecification()

    for i in range(mol.GetNumConformers()):
        conf = mol.GetConformer(i)
        single_conf_spec = pyrosetta.rosetta.std.map_std_string_numeric_xyzVector_double_t_std_allocator_std_pair_const_std_string_numeric_xyzVector_double_t()
        for idx, atm_vd in index_to_vd.items():
            rdkit_atm_pos = conf.GetAtomPosition(idx)
            single_conf_spec[ mutable_res.atom_name(atm_vd) ] = rosetta.numeric.xyzVector_double_t(rdkit_atm_pos.x, rdkit_atm_pos.y, rdkit_atm_pos.z)

        rotamers_spec.add_rotamer(single_conf_spec)

    mutable_res.rotamer_library_specification(rotamers_spec)
    return mutable_res

def mutable_res_to_res(mutable_res):
    lig_restype_non_mutable = rosetta.core.chemical.ResidueType.make( mutable_res )
    return rosetta.core.conformation.Residue( lig_restype_non_mutable, True )

def prepare_pose( smiles, pose, ref_mol ):

    # Generate RDKit Mol from SMILES
    new_mol = Chem.MolFromSmiles( smiles )

    # Determine number of different conformers based on molecule flexibility
    nrotbonds = Lipinski.NumRotatableBonds(new_mol)
    nconf = 49
    if nrotbonds > 12:
        nconf = 299
    elif nrotbonds >= 8:
        nconf = 199
    new_mol = Chem.AddHs( new_mol )

    # Align new Mol with reference Mol and generate conformers
    AllChem.ConstrainedEmbed( new_mol, ref_mol )
    AllChem.EmbedMultipleConfs(new_mol, numConfs=nconf, clearConfs = False, maxAttempts = 30)
    AllChem.AlignMolConformers(new_mol)

    # Turn RDKit Mol into Rosetta Residue
    mutres, index_to_vd, index_to_name = rdkit_to_mutable_res( new_mol )
    mutres = add_confs_to_res( new_mol, mutres, index_to_vd )
    res = mutable_res_to_res(mutres)

    # Examplified creation of 2D image with Rosetta atom names
    d = rdMolDraw2D.MolDraw2DCairo(750, 750) # or MolDraw2DSVG to get SVGs
    for i in range( new_mol.GetNumAtoms() ):
        new_mol.GetAtomWithIdx(i).SetProp('atomNote', index_to_name[i])
    AllChem.Compute2DCoords(new_mol)
    d.DrawMolecule(new_mol)
    d.FinishDrawing()
    d.WriteDrawingText('atom_map.png')

    # copy receptor pose (should not contain a ligand!)
    new_pose = pyrosetta.rosetta.core.pose.Pose()
    new_pose.detached_copy(pose)

    # Add ligand residue and update information
    new_pose.append_residue_by_jump( res, 1, "", "", True )
    new_pose.pdb_info().chain( new_pose.total_residue(), 'X' )
    new_pose.update_pose_chains_from_pdb_chains()

    return new_pose, new_mol

def add_constraints( pose, constraints ):

    for constraint in constraints:
        # Receptor atom selected via Residue index and atom name
        atom1 = pyrosetta.rosetta.core.id.AtomID( pose.residue( constraint[ 0 ] ).type().atom_index( constraint[ 1 ] ), constraint[ 0 ] )
        # ligand atom selected from last residue (aka ligand residue) and atom name
        atom2 = pyrosetta.rosetta.core.id.AtomID( pose.residue( pose.total_residue() ).type().atom_index( constraint[ 2 ] ), pose.total_residue() )
        harmonic = pyrosetta.rosetta.core.scoring.func.HarmonicFunc( constraint[3], constraint[4] )
        cst = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint( atom1, atom2, harmonic )
        pose.add_constraint(cst)

def dock( pose, all_results, mover, scfx, repeats=1 ):

    for _ in range(repeats):
        work_pose = pyrosetta.rosetta.core.pose.Pose()
        work_pose.detached_copy(pose)
        mover.apply( work_pose )

        results = {
            'pose': work_pose,
            'total_score': work_pose.energies().total_energy(),
            'atom_pair_cst': work_pose.energies().total_energies()[ rosetta.core.scoring.ScoreType.atom_pair_constraint ],
        }

        work_pose.remove_constraints()
        interface_scores = rosetta.protocols.ligand_docking.get_interface_deltas( 'X', work_pose, scfx )
        idelta_scores = {}
        score_types = [ 'if_X_fa_atr', 'if_X_fa_elec', 'if_X_fa_rep', 'if_X_fa_sol', 'if_X_hbond_bb_sc', 'if_X_hbond_sc', 'interface_delta_X', 'if_X_fa_pair' ]
        for score_type in score_types:
            idelta_scores[ score_type ] = interface_scores[ score_type ]
        

        results[ 'atom_attraction' ] = idelta_scores[ 'if_X_fa_atr' ]
        results[ 'electrostatic' ] = idelta_scores[ 'if_X_fa_elec' ]
        results[ 'atom_repulsion' ] = idelta_scores[ 'if_X_fa_rep' ]
        results[ 'solvation' ] = idelta_scores[ 'if_X_fa_sol' ]
        results[ 'hbond' ] = idelta_scores[ 'if_X_hbond_bb_sc' ] + idelta_scores[ 'if_X_hbond_sc' ]
        results[ 'delta_g' ] = idelta_scores[ 'interface_delta_X' ]
        results[ 'pairwise_energy' ] = idelta_scores[ 'if_X_fa_pair' ]

        all_results.append( results )

def analyze_results( results, rdkit_mol, max_cst_penalty ):

    best_pose = None
    best_score = 999999.99
    best_index = 0

    for idx, result in enumerate(results):
        if result[ 'delta_g' ] < best_score and result['atom_pair_cst'] < max_cst_penalty:
            best_score = result[ 'delta_g' ]
            best_pose = result[ 'pose' ]
            best_index = idx

    res_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(best_pose.total_residue())
    rmsd_calc = rosetta.core.simple_metrics.metrics.RMSDMetric()
    rmsd_calc.set_residue_selector(res_selector)
    rmsd_calc.set_comparison_pose( best_pose )

    for result in results:
        rmsd = rmsd_calc.calculate( result[ 'pose' ] )
        result[ 'rmsd' ] = rmsd

    x = []
    y = []
    for result in results:
        rmsd = result[ 'rmsd' ]
        score = result[ 'delta_g' ]
        x.append( rmsd )
        y.append( score )

    plt.scatter( x, y )
    plt.ylabel( 'Delta G (REU)' )
    plt.xlabel( 'RMSD' )
    plt.savefig( 'funnel.png' )

    best_pose.dump_pdb("constraints/best_pose.pdb")

    header = [ 'name' ]
    for key in results[0].keys():
        if key != 'pose':
            header.append( key )
    padding = 16

    weight = Descriptors.MolWt(rdkit_mol)
    hbond_acc = Descriptors.NOCount(rdkit_mol)
    hbond_don = Descriptors.NHOHCount(rdkit_mol)
    logp = Descriptors.MolLogP(rdkit_mol)
    qed = QED.qed( rdkit_mol )
    print( ''.rjust( padding ), ''.join( h.rjust( padding ) for h in ['weight', 'hbond_acc', 'hbond_don', 'logp', 'qed'] ), sep='' )
    print( 'Ligand Prop'.rjust( padding ), ''.join( f'{v:.4f}'.rjust( padding ) for v in [weight, hbond_acc, hbond_don, logp, qed] ), sep='' )
    print()

    print( "".join([h.rjust( padding ) for h in header]) )
    print( 'Best'.rjust( padding ), ''.join( f'{results[best_index][h]:.4f}'.rjust( padding ) for h in header[1:] ), sep='' )

    for i in range( len(results) ):
        print( str(i+1).rjust( padding ), ''.join( f'{results[i][h]:.4f}'.rjust( padding ) for h in header[1:] ), sep='' )

def main():

    pybel_mol = next(pybel.readfile("sdf", 'constraints/relax.sdf'))
    pybel_sdf = pybel_mol.write('sdf')

    # RDKits sanitization step automatically detects aromaticity
    # I am assuming that it is best to provide molecules in kekulized form with all hydrogens present
    ref_mol = Chem.MolFromMolBlock( pybel_sdf )

    pose = pyrosetta.pose_from_pdb( "constraints/relax_nolig.pdb" )

    smiles = 'COc5cc(OCc1ccncc1)cc6nc(c4cc(CCc2ccnnc2)c(OCCO)c(Cc3cnccn3)c4)[nH]c(=O)c56'
    lig_pose, new_mol = prepare_pose( smiles, pose, ref_mol )

    constraints = [
        [92, 'HE2', 'O5', 2.5, 0.25],
        [45, 'HH', 'N1', 2.5, 0.25],
        [87, 'OH', 'O4', 2.5, 0.25],
    ]

    add_constraints( lig_pose, constraints )

    protocol_xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file("constraints/transform_std.xml")
    mover = protocol_xml.get_mover("ParsedProtocol")
    scfx = protocol_xml.get_score_function("hard_rep")

    results = []
    dock( lig_pose, results, mover, scfx, repeats=150 )

    analyze_results( results, new_mol, max_cst_penalty=10 )


if __name__ == "__main__":
    main()