
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.dirname("DiffDock/")))

from DiffDock.utils.sampling import randomize_position, sampling
from DiffDock.utils.utils import get_model
from argparse import Namespace
from functools import partial
from DiffDock.utils.diffusion_utils import t_to_sigma as t_to_sigma_compl, get_t_schedule
from DiffDock.utils.inference_utils import InferenceDataset, set_nones
from torch_geometric.loader import DataLoader
from rdkit.Chem import RemoveAllHs
from rdkit.Geometry import Point3D

import numpy as np

import torch
import yaml
import copy
import pandas as pd

import warnings


warnings.filterwarnings("ignore", category=UserWarning,
                        message="The TorchScript type system doesn't support instance-level annotations on empty non-base types in `__init__`")

inference_steps = 20
model_dir = 'DiffDock/workdir/v1.1/score_model'
conf_model_dir = 'DiffDock/workdir/v1.1/confidence_model'
out_dir = 'output'
samples_per_complex = 5
center = np.array([140.2520, 126.5740, 133.8730])

with open(model_dir+'/model_parameters.yml') as f:
    score_model_args = Namespace(**yaml.full_load(f))

with open(conf_model_dir + '/model_parameters.yml') as f:
    confidence_args = Namespace(**yaml.full_load(f))

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print("----DiffDock is running on", device)

t_to_sigma = partial(t_to_sigma_compl, args=score_model_args)
tr_schedule = get_t_schedule(inference_steps=inference_steps, sigma_schedule='expbeta')

model = get_model(score_model_args, device, t_to_sigma=t_to_sigma, no_parallel=True, old=False)
state_dict = torch.load(model_dir + '/best_ema_inference_epoch_model.pt', map_location=torch.device('cpu'))
model.load_state_dict(state_dict, strict=True)
model = model.to(device)
model.eval()

confidence_model = get_model(confidence_args, device, t_to_sigma=t_to_sigma, no_parallel=True,
                                    confidence_mode=True, old=True)
state_dict = torch.load(conf_model_dir + '/best_model_epoch75.pt', map_location=torch.device('cpu'))
confidence_model.load_state_dict(state_dict, strict=True)
confidence_model = confidence_model.to(device)
confidence_model.eval()

def diffuse_ligand(smiles):

    complex_name_list = set_nones(['8ef6'])
    protein_path_list = set_nones(['input/8ef6_0001_R.pdb'])
    ligand_description_list = set_nones([smiles])
    protein_sequence_list = set_nones([None])

    test_dataset = InferenceDataset(out_dir=out_dir, complex_names=complex_name_list, protein_files=protein_path_list,
                                    ligand_descriptions=ligand_description_list, protein_sequences=protein_sequence_list,
                                    lm_embeddings=True,
                                    receptor_radius=score_model_args.receptor_radius, remove_hs=score_model_args.remove_hs,
                                    c_alpha_max_neighbors=score_model_args.c_alpha_max_neighbors,
                                    all_atoms=score_model_args.all_atoms, atom_radius=score_model_args.atom_radius,
                                    atom_max_neighbors=score_model_args.atom_max_neighbors,
                                    knn_only_graph=False if not hasattr(score_model_args, 'not_knn_only_graph') else not score_model_args.not_knn_only_graph)
    test_loader = DataLoader(dataset=test_dataset, batch_size=1, shuffle=False)

    confidence_test_dataset = \
            InferenceDataset(out_dir=out_dir, complex_names=complex_name_list, protein_files=protein_path_list,
                             ligand_descriptions=ligand_description_list, protein_sequences=protein_sequence_list,
                             lm_embeddings=True,
                             receptor_radius=confidence_args.receptor_radius, remove_hs=confidence_args.remove_hs,
                             c_alpha_max_neighbors=confidence_args.c_alpha_max_neighbors,
                             all_atoms=confidence_args.all_atoms, atom_radius=confidence_args.atom_radius,
                             atom_max_neighbors=confidence_args.atom_max_neighbors,
                             precomputed_lm_embeddings=test_dataset.lm_embeddings,
                             knn_only_graph=False if not hasattr(score_model_args, 'not_knn_only_graph') else not score_model_args.not_knn_only_graph)
    
    for idx, orig_complex_graph in enumerate(test_loader):
        confidence_complex_graph = confidence_test_dataset[idx]
        confidence_data_list = [copy.deepcopy(confidence_complex_graph) for _ in range(samples_per_complex)]
        data_list = [copy.deepcopy(orig_complex_graph) for _ in range(samples_per_complex)]

        randomize_position(data_list, score_model_args.no_torsion, False, score_model_args.tr_sigma_max,
                               initial_noise_std_proportion=1.4601642460337794,
                               choose_residue=False)

        lig = orig_complex_graph.mol[0]

        # run reverse diffusion
        print("Start inference")
        n_fails = 0
        while n_fails < 5:
            try:
                data_list, confidence = sampling(data_list=data_list, model=model,
                                                    inference_steps=inference_steps-1,
                                                    tr_schedule=tr_schedule, rot_schedule=tr_schedule, tor_schedule=tr_schedule,
                                                    device=device, t_to_sigma=t_to_sigma, model_args=score_model_args,
                                                    visualization_list=None, confidence_model=confidence_model,
                                                    confidence_data_list=confidence_data_list, confidence_model_args=confidence_args,
                                                    batch_size=1, no_final_step_noise=True,
                                                    temp_sampling=[1.170050527854316, 2.06391612594481,
                                                                7.044261621607846],
                                                    temp_psi=[0.727287304570729, 0.9022615585677628, 0.5946212391366862],
                                                    temp_sigma_data=[0.9299802531572672, 0.7464326999906034,
                                                                    0.6943254174849822])
                break
            except:
                n_fails += 1
                print("Inference failed", n_fails, "times")
        print("Done with inference")
        
        ligand_pos = np.asarray([complex_graph['ligand'].pos.cpu().numpy() + orig_complex_graph.original_center.cpu().numpy() for complex_graph in data_list])

        # reorder predictions based on confidence output
        if confidence is not None and isinstance(confidence_args.rmsd_classification_cutoff, list):
            confidence = confidence[:, 0]
        if confidence is not None:
            confidence = confidence.cpu().numpy()
            re_order = np.argsort(confidence)[::-1]
            confidence = confidence[re_order]
            ligand_pos = ligand_pos[re_order]

        distances = []
        for positions in ligand_pos:
            pos = np.sum(positions, 0) / len(positions)
            distances.append(np.linalg.norm(pos-center))

        best_result_idx = len(ligand_pos) + 1
        for i in range(len(ligand_pos)):
            if distances[i] < 10.0:
                best_result_idx = i
                print("Picked result", best_result_idx, "with distance", distances[i], "and confidence", confidence[i])
                break

        if best_result_idx == len(ligand_pos) + 1:
            best_result_idx = 0
            print("Failed to dock into pocket. Continue with distance", distances[best_result_idx], "and best confidence", confidence[best_result_idx])

        mol_pred = copy.deepcopy(lig)
        if score_model_args.remove_hs: mol_pred = RemoveAllHs(mol_pred)
        conf = mol_pred.GetConformer()

        for i in range(mol_pred.GetNumAtoms()):
            x,y,z = ligand_pos[best_result_idx].astype(np.double)[i]
            conf.SetAtomPosition(i,Point3D(x,y,z))

    return mol_pred, distances[best_result_idx]