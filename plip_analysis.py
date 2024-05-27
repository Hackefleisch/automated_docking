from plip.structure.preparation import PDBComplex
from plip.exchange.report import StructureReport
from plip.basic.remote import VisualizerData

from plip.visualization.pymol import PyMOLVisualizer
from pymol import cmd
import pymol
import sys

def initialize_pymol(options):
    """Initializes PyMOL"""
    # Pass standard arguments of function to prevent PyMOL from printing out PDB headers (workaround)
    pymol.finish_launching(args=['pymol', options, '-K'])
    pymol.cmd.reinitialize()

def start_pymol(quiet=False, options='-p', run=False):
    """Starts up PyMOL and sets general options. Quiet mode suppresses all PyMOL output.
    Command line options can be passed as the second argument."""
    pymol.pymol_argv = ['pymol', '%s' % options] + sys.argv[1:]
    if run:
        initialize_pymol(options)
    if quiet:
        pymol.cmd.feedback('disable', 'all', 'everything')

def visualize_in_pymol(plcomplex):
    """Visualizes the given Protein-Ligand complex at one site in PyMOL."""

    vis = PyMOLVisualizer(plcomplex)

    #####################
    # Set everything up #
    #####################

    pdbid = plcomplex.pdbid
    lig_members = plcomplex.lig_members
    chain = plcomplex.chain

    ligname = vis.ligname
    hetid = plcomplex.hetid

    metal_ids = plcomplex.metal_ids
    metal_ids_str = '+'.join([str(i) for i in metal_ids])

    ########################
    # Basic visualizations #
    ########################

    start_pymol(run=True, options='-pcq', quiet=False)
    vis.set_initial_representations()

    cmd.load(plcomplex.sourcefile)
    cmd.frame(1)
    current_name = cmd.get_object_list(selection='(all)')[0]

    cmd.set_name(current_name, pdbid)
    cmd.hide('everything', 'all')
    cmd.select(ligname, 'resn %s and chain %s and resi %s*' % (hetid, chain, plcomplex.position))

    # Visualize and color metal ions if there are any
    if not len(metal_ids) == 0:
        vis.select_by_ids(ligname, metal_ids, selection_exists=True)
        cmd.show('spheres', 'id %s and %s' % (metal_ids_str, pdbid))

    # Additionally, select all members of composite ligands
    if len(lig_members) > 1:
        for member in lig_members:
            resid, chain, resnr = member[0], member[1], str(member[2])
            cmd.select(ligname, '%s or (resn %s and chain %s and resi %s)' % (ligname, resid, chain, resnr))

    cmd.show('sticks', ligname)
    cmd.color('myblue')
    cmd.color('myorange', ligname)
    cmd.util.cnc('all')
    if not len(metal_ids) == 0:
        cmd.color('hotpink', 'id %s' % metal_ids_str)
        cmd.hide('sticks', 'id %s' % metal_ids_str)
        cmd.set('sphere_scale', 0.3, ligname)
    cmd.deselect()

    vis.make_initial_selections()

    vis.show_hydrophobic()  # Hydrophobic Contacts
    vis.show_hbonds()  # Hydrogen Bonds
    vis.show_halogen()  # Halogen Bonds
    vis.show_stacking()  # pi-Stacking Interactions
    vis.show_cationpi()  # pi-Cation Interactions
    vis.show_sbridges()  # Salt Bridges
    vis.show_wbridges()  # Water Bridges
    vis.show_metal()  # Metal Coordination

    vis.refinements()

    vis.zoom_to_ligand()

    vis.selections_group()
    vis.additional_cleanup()
    
    vis.save_picture("output/", "picture")
    cmd.mstop()

    vis.selections_cleanup()
    cmd.load(plcomplex.sourcefile)
    # Create a selection for atoms within 5 angstroms of the ligand
    cmd.select('near_ligand', 'byres (resn UNK around 5)')

    # Show surface for the selection
    cmd.show('surface', 'near_ligand')
    cmd.show('sticks', 'near_ligand')
    cmd.hide('cartoon', 'near_ligand')
    cmd.set('transparency', 0.5, 'near_ligand')

    cmd.select('far_ligand', 'not near_ligand and not resn UNK')
    cmd.hide('everything', 'far_ligand')
    cmd.deselect()
    
    cmd.color('myblue')
    cmd.color('myorange', "resn UNK")
    cmd.util.cnc('all')

    cmd.center("resn UNK")
    cmd.orient("resn UNK")
    cmd.zoom("resn UNK", 3)
    cmd.origin("resn UNK")

    cmd.save("output/pymol_session.pse")

def analyse_binding(mol):
    site = list(mol.interaction_sets.values())[0]
    
    n_hydroph = len(site.hydrophobic_contacts)
    n_hbond = len(site.hbonds_pdon + site.hbonds_ldon)
    n_saltbridge = len(site.saltbridge_lneg + site.saltbridge_pneg)
    n_pistack = len(site.pistacking)
    n_pication = len(site.pication_laro + site.pication_paro)
    n_halogen = len(site.halogen_bonds)

    # check if res 84 participates in hbond or saltbridge
    interact_84 = False
    for sb in site.saltbridge_lneg + site.saltbridge_pneg:
        if sb.resnr == 84:
            interact_84 = True
            break

    #if not interact_84:
    #    for hbond in site.hbonds_pdon + site.hbonds_ldon:
    #        if hbond.resnr == 84:
    #            interact_84 = True
    #            break

    return n_hydroph, n_hbond, n_saltbridge, n_pistack, n_pication, n_halogen, interact_84

def run_analysis():
    mol = PDBComplex()
    mol.load_pdb('output/tmp.pdb')
    for ligand in mol.ligands:
        mol.characterize_complex(ligand)

    streport = StructureReport(mol, outputprefix="")
    with open("output/tmp_report.txt", 'w') as file:
        file.write("\n".join(streport.txtreport))

    complexes = [VisualizerData(mol, site) for site in sorted(mol.interaction_sets)
                    if not len(mol.interaction_sets[site].interacting_res) == 0]
    complexes[0].sourcefile = 'output/tmp.pdb'
    print(complexes[0].sourcefile)
    [visualize_in_pymol(plcomplex) for plcomplex in complexes]

    return analyse_binding(mol)

def main():
    run_analysis()

if __name__ == "__main__":
    main()