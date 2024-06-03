import panel as pn
from panel_chemistry.widgets import JSMEEditor
from bokeh.models.widgets.tables import NumberFormatter
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import pandas as pd
import os
import shutil
import torch


import ligand_docking
import plip_analysis

pn.extension("jsme", sizing_mode="stretch_width")

editor = JSMEEditor(value="CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O", height=500, format="smiles")

#editor.servable()

#pn.extension('tabulator', css_files=["https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.2/css/all.min.css"])
pn.extension('tabulator')

pose = ligand_docking.load_pose( "input/8ef6_0001_R.pdb" )
protocol, scfx = ligand_docking.create_protocol("input/transform_repack.xml", ["high_res_docker", "final"], ["hard_rep"])
complex = None

cwd = os.getcwd()
print("CWD:", cwd)
ligand_image = 'output/tmp_2dmol.png'

# current data as a Pandas DataFrame
current_data = {
    'Name': [ None ], 
    'Weight': [ None ], 
    'Hbond acc': [ None ], 
    'Hbond don': [ None ], 
    'LogP': [ None ], 
    'Score': [ None ],
    'Pocket': [ None ],
    '84': [ None ],
    'Hydrophobic': [ None ],
    'Hbond': [ None ],
    'Saltbridge': [ None ],
    'Pistack': [ None ],
    'Pication': [ None ],
    'Halogen': [ None ],
    'Smiles' : [ None ],
}
current_df = pd.DataFrame(current_data)
current_df.set_index('Name', inplace=True)

saved_data = {
    'Name': [], 
    'Weight': [  ], 
    'Hbond acc': [  ], 
    'Hbond don': [  ], 
    'LogP': [  ], 
    'Score': [  ],
    'Pocket': [  ],
    '84': [  ],
    'Hydrophobic': [  ],
    'Hbond': [  ],
    'Saltbridge': [  ],
    'Pistack': [  ],
    'Pication': [  ],
    'Halogen': [  ],
    'image': [  ],
    'Smiles' : [  ],
}
if os.path.exists("results/docking.save"):
    saved_df = pd.read_pickle("results/docking.save")
else:
    saved_df = pd.DataFrame(saved_data)
    saved_df.set_index('Name', inplace=True)

# Create a Tabulator table
image_formatter = {
    'type' : 'image',
    'height' : "50px",
    'width' : "auto",
    'urlPrefix' : "",
    'urlSuffix' : "",
}
text_formatter = {
    'type' : 'text'
}

def color_red_5(val):
    color = 'red' if val == None or val > 5 else 'green'
    return 'color: %s' % color

def color_red_10(val):
    color = 'red' if val == None or val > 10 else 'green'
    return 'color: %s' % color

def color_red_500(val):
    color = 'red' if val == None or val > 500 else 'green'
    return 'color: %s' % color

formatter = {
    'image' : image_formatter,
    'Name': text_formatter, 
    'Pocket': {'type': 'tickCross'}, 
    '84': {'type': 'tickCross'}, 
    'Score':  NumberFormatter(format='0.0000')
}

table_buttons={
    'Load': '<img src="https://cdn4.iconfinder.com/data/icons/basic-user-interface-elements/700/export-share-upload-512.png" width="20" height="20"></img>',
    '_Delete': '<img src="https://cdn3.iconfinder.com/data/icons/font-awesome-regular-1/512/trash-can-512.png" width="20" height="20"></img>'
}

current_mol_table = pn.widgets.Tabulator(current_df, formatters=formatter, embed_content=False, layout='fit_columns', disabled=True, hidden_columns =['Smiles'])
current_mol_table.sizing_mode = 'stretch_width'
current_mol_table.style.applymap(color_red_5, subset=['Hbond don', 'LogP'])
current_mol_table.style.applymap(color_red_10, subset=['Hbond acc'])
current_mol_table.style.applymap(color_red_500, subset=['Weight'])

saved_mol_table = pn.widgets.Tabulator(saved_df, formatters=formatter, embed_content=False, pagination='remote', page_size=30, layout='fit_columns', buttons=table_buttons, disabled=True, hidden_columns =['Smiles'])
saved_mol_table.sizing_mode = 'stretch_width'
saved_mol_table.style.applymap(color_red_5, subset=['Hbond don', 'LogP'])
saved_mol_table.style.applymap(color_red_10, subset=['Hbond acc'])
saved_mol_table.style.applymap(color_red_500, subset=['Weight'])


def table_click(event):
    if event.column == "_Delete":
        name = saved_df.index.to_list()[event.row]
        print("Delete", name)
        saved_df.drop(name, inplace=True)
        print(saved_df)
        saved_mol_table.value = saved_df
        saved_df.to_pickle("results/docking.save")
        shutil.rmtree("results/" + name)
    elif event.column == "Load":
        name = saved_df.index.to_list()[event.row]
        print("Load", name)
        load_current_data(name)

saved_mol_table.on_click(table_click)

mol2d_image = pn.pane.PNG(None, height=150)
interaction_image = pn.pane.PNG(None, height=350)

ligand_name_input = pn.widgets.TextInput(value='', placeholder='Ligand name')

def analyze_mol(event):
    global editor, current_data, current_df, current_mol_table
    if current_data[ 'Weight' ] != [ None ]:
        info_field.object = "Reset values befor new calculation"
        return
    smiles = editor.value
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        AllChem.Compute2DCoords(mol)
        Draw.MolToFile(mol, ligand_image)
        current_data[ 'Weight' ] = [ Descriptors.MolWt(mol) ]
        current_data[ 'Hbond acc' ] = [ Descriptors.NOCount(mol) ]
        current_data[ 'Hbond don' ] = [ Descriptors.NHOHCount(mol) ]
        current_data[ 'LogP' ] = [ Descriptors.MolLogP(mol) ] 
        current_data[ 'Smiles' ] = smiles
        current_df = pd.DataFrame(current_data)
        current_df.set_index('Name', inplace=True)
        current_mol_table.value = current_df
        mol2d_image.object = None
        mol2d_image.object = ligand_image

def dock_ligand(event):
    global pose, protocol, complex, scfx, editor, current_data, current_df, current_mol_table
    smiles = editor.value
    if current_data[ 'Score' ] != [ None ]:
        info_field.object = "Reset values befor new calculation"
        return
    info_field.object = "Running calculation..."
    docker = ligand_docking.FullDocker(smiles, pose, protocol, scfx)
    finished = False
    while not finished:
        info_field.object = "Running calculation... " + docker.next_step
        finished = docker.run()
    complex = docker.best_complex
    score = docker.best_score
    distance = docker.distance
    info_field.object = "Running calculation... Done. Ligand has distance of " + str(distance)
    current_data[ 'Score' ] = [ score ] 
    current_data[ 'Pocket' ] = [ distance < 8.0 ] 
    current_df = pd.DataFrame(current_data)
    current_df.set_index('Name', inplace=True)
    current_mol_table.value = current_df
    complex.dump_pdb("output/tmp.pdb")

def analyze_pose(event):
    if current_data[ 'Hydrophobic' ] != [ None ]:
        info_field.object = "Reset values befor new calculation"
        return
    n_hydroph, n_hbond, n_saltbridge, n_pistack, n_pication, n_halogen, interact_84 = plip_analysis.run_analysis()
    current_data[ '84' ] = [ interact_84 ] 
    current_data[ 'Hydrophobic' ] = [ n_hydroph ] 
    current_data[ 'Hbond' ] = [ n_hbond ] 
    current_data[ 'Saltbridge' ] = [ n_saltbridge ] 
    current_data[ 'Pistack' ] = [ n_pistack ] 
    current_data[ 'Pication' ] = [ n_pication ] 
    current_data[ 'Halogen' ] = [ n_halogen ] 
    current_df = pd.DataFrame(current_data)
    current_df.set_index('Name', inplace=True)
    current_mol_table.value = current_df
    interaction_image.object = None
    interaction_image.object = 'output/picture.png'

def save_results(event):
    global saved_df, current_data, current_df, current_mol_table
    lig_name = ligand_name_input.value
    if lig_name == "":
        info_field.object = "Name is empty."
        return
    if lig_name in saved_df.index:
        info_field.object = "Name alreay in use. Select a differen one."
        return
    current_data[ 'Name' ] = [ lig_name ]
    current_df = pd.DataFrame(current_data)
    current_df.set_index('Name', inplace=True)
    current_mol_table.value = current_df
    if not current_df.isnull().values.any():
        saved_df = pd.concat([current_df, saved_df])
        saved_df[ 'image' ][ lig_name ] = [ "http://localhost:8000/results/" + lig_name + "/mol2d.png" ]
        saved_mol_table.value = saved_df
        reset_current_data()
        move_files(lig_name)
        saved_df.to_pickle("results/docking.save")
        ligand_name_input.value = ""
    else:
        info_field.object = "Finish calculations before saving"

def move_files(name):
    os.makedirs("results/" + name, exist_ok=True)

    source_path = 'output/picture.png'
    destination_path = 'results/' + name + '/interactions.png'
    shutil.copy(source_path, destination_path)

    source_path = 'output/pymol_session.pse'
    destination_path = 'results/' + name + '/pymol_session.pse'
    shutil.copy(source_path, destination_path)

    source_path = 'output/tmp_2dmol.png'
    destination_path = 'results/' + name + '/mol2d.png'
    shutil.copy(source_path, destination_path)

    source_path = 'output/tmp_report.txt'
    destination_path = 'results/' + name + '/report.txt'
    shutil.copy(source_path, destination_path)

    source_path = 'output/tmp.pdb'
    destination_path = 'results/' + name + '/' + name + '.pdb'
    shutil.copy(source_path, destination_path)

def reset_current_data(event=None):
    global current_data, current_df
    current_data = {
        'Name': [ None ], 
        'Weight': [ None ], 
        'Hbond acc': [ None ], 
        'Hbond don': [ None ], 
        'LogP': [ None ], 
        'Score': [ None ],
        'Pocket': [ None ],
        '84': [ None ],
        'Hydrophobic': [ None ],
        'Hbond': [ None ],
        'Saltbridge': [ None ],
        'Pistack': [ None ],
        'Pication': [ None ],
        'Halogen': [ None ],
        'Smiles': [ None ],
    }
    current_df = pd.DataFrame(current_data)
    current_df.set_index('Name', inplace=True)
    current_mol_table.value = current_df
    interaction_image.object = None
    mol2d_image.object = None

def load_current_data(mol_name):
    global current_df, saved_df
    if current_df.isnull().all().all():
        current_name = current_df.index.to_list()[0]
        current_df.loc[current_name] = saved_df.loc[mol_name]
        #current_df.set_index('Name', inplace=True)
        current_mol_table.value = current_df
        interaction_image.object = 'results/' + mol_name + '/interactions.png'
        mol2d_image.object = 'results/' + mol_name + '/mol2d.png'
        editor.value = saved_df.loc[mol_name][ 'Smiles' ]
    else:
        info_field.object = "Clear results before loading"

# Create a button widget
draw_button = pn.widgets.Button(name='Generate Image', button_type='primary')
dock_button = pn.widgets.Button(name='Dock Ligand', button_type='primary')
analyze_button = pn.widgets.Button(name='Analyze Pose', button_type='primary')
save_button = pn.widgets.Button(name='Save results', button_type='primary')
reset_button = pn.widgets.Button(name='Reset results', button_type='primary')

# Link the button to the function
draw_button.on_click(analyze_mol)
dock_button.on_click(dock_ligand)
analyze_button.on_click(analyze_pose)
save_button.on_click(save_results)
reset_button.on_click(reset_current_data)

info_field = pn.pane.Markdown()

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
cuda_field = pn.pane.Markdown("Torch device: " + str(device))

# Create and display the Panel layout
layout = pn.Column(
    pn.Row(editor, pn.Column(mol2d_image, interaction_image, align='center')), 
    current_mol_table, 
    info_field,
    pn.Row(reset_button, draw_button, dock_button, analyze_button, ligand_name_input, save_button), 
    saved_mol_table,
    cuda_field,
)
layout.servable()