import panel as pn
from panel_chemistry.widgets import JSMEEditor
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import PIL
from io import BytesIO

import ligand_docking
import plip_analysis

pn.extension("jsme", sizing_mode="stretch_width")

editor = JSMEEditor(value="CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O", height=500, format="smiles")

editor.servable()

pose = ligand_docking.load_pose( "input/8ef6_0001_R.pdb" )
protocol, scfx = ligand_docking.create_protocol("input/transform_repack.xml", ["high_res_docker", "final"], ["hard_rep"])
complex = None

def draw_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)  # You can adjust size here
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        byte_string = drawer.GetDrawingText()
        return PIL.Image.open(BytesIO(byte_string))
    else:
        return None
    
def update_image(event):
    smiles = editor.value
    image_pane.object = draw_mol(smiles)

def dock_ligand(event):
    global pose, protocol, complex, scfx, editor
    smiles = editor.value
    complex, score = ligand_docking.full_docking(smiles, pose, protocol, scfx)
    print(complex.scores)
    print(score)
    complex.dump_pdb("output/tmp.pdb")

def analyze_pose(event):
    plip_analysis.run_analysis()

image_pane = pn.pane.PNG()

# Create a button widget
draw_button = pn.widgets.Button(name='Generate Image', button_type='primary')
dock_button = pn.widgets.Button(name='Dock Ligand', button_type='primary')
analyze_button = pn.widgets.Button(name='Analyze Pose', button_type='primary')

# Link the button to the function
draw_button.on_click(update_image)
dock_button.on_click(dock_ligand)
analyze_button.on_click(analyze_pose)

# Create and display the Panel layout
layout = pn.Column(editor, pn.Row(draw_button, dock_button, analyze_button), image_pane)
layout.show()