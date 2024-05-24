import panel as pn
from panel_chemistry.widgets import JSMEEditor
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import PIL
from io import BytesIO

import ligand_docking


pn.extension("jsme", sizing_mode="stretch_width")

editor = JSMEEditor(value="CC", height=500, format="smiles")

editor.servable()

pose = None
rosetta_protocol = ligand_docking.create_protocol("input/transform_repack.xml")

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
    global pose
    smiles = editor.value
    res = ligand_docking.smiles_to_res(smiles)
    pose = ligand_docking.load_pose("input/protein_target.pdb")
    pose = ligand_docking.create_complex(pose, res)
    rosetta_protocol.apply(pose)
    print(pose.scores)

def save_pose(event):
    global pose
    pose.dump_pdb("test.pdb")

image_pane = pn.pane.PNG()

# Create a button widget
draw_button = pn.widgets.Button(name='Generate Image', button_type='primary')
dock_button = pn.widgets.Button(name='Dock Ligand', button_type='primary')
save_button = pn.widgets.Button(name='Save Pose', button_type='primary')

# Link the button to the function
draw_button.on_click(update_image)
dock_button.on_click(dock_ligand)
save_button.on_click(save_pose)

# Create and display the Panel layout
layout = pn.Column(editor, pn.Row(draw_button, dock_button, save_button), image_pane)
layout.show()