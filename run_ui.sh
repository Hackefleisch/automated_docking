#!/bin/bash

trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

git pull

python -m http.server 8000  &

eval "$(conda shell.bash hook)"
conda activate p6docking38
panel serve panel.py 