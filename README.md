# DREAMM
Predicting protein-membrane interfaces using ensemple machine learning

## Installation
conda create -n dreamm moleculekit htmd-pdb2pqr biopython prody mdanalysis scikit-learn dssp msms hhsuite lightgbm mlxtend -c acellera -c anaconda -c insilichem -c conda-forge -c salilab -c bioconda

conda activate dreamm
pip install freesasa

chmod +x dreamm/dreamm/dreamm.py
