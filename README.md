# DREAMM
Predicting protein-membrane interfaces using ensemple machine learning


## Installation
conda create -n dreamm moleculekit htmd-pdb2pqr biopython prody mdanalysis scikit-learn dssp msms hhsuite lightgbm mlxtend -c acellera -c anaconda -c insilichem -c conda-forge -c salilab -c bioconda

conda activate dreamm

pip install freesasa

git clone https://github.com/zoecournia/DREAMM

chmod +x DREAMM/dreamm/dreamm.py

export DREAMM="/home/USER/DREAMM/dreamm" #Replace USER with your username

alias dreamm="$DREAMM/dreamm.py"


## Usage
dreamm -h

Input Settings:
  -i INPUT, --input INPUT
                        The 4 letter PDB code or the path to the PDB file
  -c CHAIN [CHAIN ...], --chain CHAIN [CHAIN ...]
                        The chain or chains. Default is None
  -d DATABASE, --database DATABASE
                        The path to the Uniclust30_2018_8 database.
  -p PROCESSES, --processes PROCESSES
                        The number of CPUs. Default is All minus 1
                        
## Example
dreamm -i 5BZZ -c A -d /media/mydata/Databases/uniclust30_2018_08_hhsuite/uniclust30_2018_08 -p 16

