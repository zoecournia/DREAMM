# About DREAMM
DREAMM provides a fast and robust prediction of protein-membrane interfaces for peripheral membrane proteins using ensemble machine learning.

## Getting Started
### Prerequisites
We recommend installing Miniconda on your machine to better manage python packages and environments.

We recommend installing DREAMM in a new conda environment (see <a href="#Installation">Installation</a>). <a href="http://example.com" target="_blank">http://example.com</a>

Java is also necessary to be installed on your machine as DREAMM utilizes also [ProtDCal](https://protdcal.zmb.uni-due.de/){:target="_blank"}, which is written in java, for feature extraction.

Lastly, you have to download [Uniclust30_2018_08_hhsuite](http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/){:target="_blank" rel="noopener"} database, which is necessary to generate conservation scores (may take a while to download).

### Installation
1. Create a new conda environment, installing all necessary python libraries
```
conda create -n dreamm moleculekit htmd-pdb2pqr biopython prody mdanalysis scikit-learn dssp msms hhsuite lightgbm mlxtend -c acellera -c anaconda -c insilichem -c conda-forge -c salilab -c bioconda
```
2. Activate the environment
```
conda activate dreamm
```
3. Install freesasa using ```pip```
```
pip install freesasa
```
4. Clone the DREAMM repo
```
git clone https://github.com/zoecournia/DREAMM
```
5. Make dreamm.py executable
```
chmod +x DREAMM/dreamm/dreamm.py
```
6. Export and alias DREAMM
```
export DREAMM="/home/USER/DREAMM/dreamm" #Replace USER with your username 
alias dreamm="$DREAMM/dreamm.py"
```

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

