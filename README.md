# About DREAMM
DREAMM provides a fast and robust prediction of protein-membrane interfaces for peripheral membrane proteins using ensemble machine learning. DREAMM is also available as a [web-server](https://dreamm.ni4os.eu/).

## Getting Started
### Prerequisites
We recommend installing Miniconda on your machine to better manage python packages and environments.

We recommend installing DREAMM in a new conda environment (see <a href="#Installation">Installation</a>). 

Java is also necessary to be installed on your machine as DREAMM utilizes also [ProtDCal](https://protdcal.zmb.uni-due.de/), which is written in java, for feature extraction.

Lastly, you have to download [Uniclust30_2018_08_hhsuite](http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/) database, which is necessary to generate conservation scores (may take a while to download).

### Installation
1. Create a new conda environment, installing all necessary python libraries
```
(base) user@computer:~$ conda create -n dreamm moleculekit pdb2pqr biopython prody mdanalysis scikit-learn dssp msms hhsuite lightgbm mlxtend -c acellera -c anaconda -c insilichem -c conda-forge -c salilab -c bioconda
```
2. Activate the environment
```
(base) user@computer:~$ conda activate dreamm
```
3. Install freesasa using ```pip```
```
(dreamm) user@computer:~$ pip install freesasa
```
4. Clone the DREAMM repo
```
(dreamm) user@computer:~$ git clone https://github.com/zoecournia/DREAMM
```
5. Make dreamm.py executable
```
(dreamm) user@computer:~$ chmod +x DREAMM/dreamm/dreamm.py
```
6. To use DREAMM from any directory, place the following in the end of your .bashrc file
```
export DREAMM="/home/USER/DREAMM/dreamm" #Replace USER with your username 
alias dreamm="$DREAMM/dreamm.py"
```

## Usage
```
(dreamm) user@computer:~$ dreamm -h
```

Input Parameters:
```
  -i INPUT, --input INPUT                           The 4 letter PDB code or the path to the PDB file
  -c CHAIN [CHAIN ...], --chain CHAIN [CHAIN ...]   The chain or chains. Default is None
  -d DATABASE, --database DATABASE                  The path to the Uniclust30_2018_8 database.
  -p PROCESSES, --processes PROCESSES               The number of CPUs. Default is All minus 1
 ```

## Example
```
(dreamm) user@computer:~$ dreamm -i 5BZZ -c A -d /media/mydata/Databases/uniclust30_2018_08_hhsuite/uniclust30_2018_08/ -p 16
```
Results are output to the terminal and written to DREAMM/dreamm/outputs/prepared/fixed/5BZZ.csv, e.g.,
```
The residues:
 A  L   42  0
 A  K  263  0
 A  M  264  0
are predicted to insert the membrane.
```
where the first column indicates the chain, the second the one-letter code amino acid type, the third the resid, and the fourth indicates if the prediction is near N- or C-terminal, or near missing loops, which might be a possible false prediction (0 for no, and 1 for yes).

## License
Distributed under the GPL-3.0 License. See `LICENSE` for more information.


## Issues
For any bugs or questions on usage feel free to use the issue tracker of this github repo.

## Citing DREAMM

If you use DREAMM in your publication please cite:

Alexios Chatzigoulas and Zoe Cournia.
**Predicting protein–membrane interfaces of peripheral membrane proteins using ensemble machine learning.**
*Briefings in Bioinformatics*, **2022**, bbab518.
[doi:10.1093/bib/bbab518](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbab518/6527274)
