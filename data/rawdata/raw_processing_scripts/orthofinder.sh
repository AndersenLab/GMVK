#!/bin/bash
# $ conda create -n orthofinder_env -y -c bioconda orthofinder

# Activate env with OrthoFinder
source activate orthofinder_env

#run orthofinder
orthofinder -f data/rawdata/ortho_prot/ -og -t 12