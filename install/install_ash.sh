#!/bin/bash
"""
5/15/24

Ash script reinstallation to fix issue 
where packages are not being found by different modules

here, conda is mamba
"""

conda create -n ash_reinstall python=3.11 git pip ipython
conda activate ash_reinstall
pip install git+https://github.com/RagnarB83/ash.git@NEW
conda install -c conda-forge openmm pdbfixer

