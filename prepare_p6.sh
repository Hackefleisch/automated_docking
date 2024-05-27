#!/bin/bash

# init

# uncomment if you want to to use a local version pyrosetta
#cp pyrosetta* /opt/conda/conda/pkgs/ -r

mkdir ~/Desktop/P6
cd ~/Desktop/P6
git clone https://github.com/Hackefleisch/automated_docking.git
cd automated_docking
git clone https://github.com/gcorso/DiffDock.git

conda env create -f environment.yml 

python -m http.server 8000  &

# run
eval "$(conda shell.bash hook)"
conda activate p6docking38
panel serve panel.py 
