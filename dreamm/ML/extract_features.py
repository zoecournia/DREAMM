# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 11:24:27 2019

@author: alexis
"""

import os, sys, subprocess, urllib.request, warnings
from multiprocessing import Pool
from Bio.PDB import *
import pandas as pd
import freesasa
from prody.dynamics import *
from prody.proteins import *
import numpy as np
from moleculekit.tools.hhblitsprofile import getSequenceProfile
from moleculekit.molecule import Molecule
from moleculekit.tools.preparation import proteinPrepare
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
warnings.filterwarnings("ignore")


substitutions = {
        'AR0':'ARG', 'ASH':'ASP', 'CYM':'CYS', 'CYX':'CYS', 'GLH':'GLU', 'HIE':'HIS', 'HSE':'HIS', 'HID':'HIS', 'HSD':'HIS', 'HIP':'HIS', 'DSG':'ASN',
        'HSP':'HIS', 'LYN':'LYS', 'TYM':'TYR', '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'TRO':'TRP',
        'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG', 'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'TRG':'LYS',
        'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA', 'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'TPQ':'ALA',
        'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS', 'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'TPO':'THR',
        'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS', 'MEN':'ASN', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR', 'TPL':'TRP',
        'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA', 'TIH':'ALA',
        'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP', 'SVA':'SER',
        'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY', 'STY':'TYR',
        'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO', 'SOC':'CYS',
        'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'SMC':'CYS',
        'MHS':'HIS', 'MIS':'SER', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU', 'SHR':'LYS',
        'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE', 'SHC':'CYS',
        'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS', 'SET':'SER',
        'SEL':'SER', 'SEP':'SER', 'DSG':'ASN', 'MED':'MET', '2TL':'THR', 'CYSF':'CYS', 'MLY':'LYS', 'CAF':'CYS', 'YCM':'CYS', 'ALS':'ALA', 'PPN':'PHE',
        'SEC':'CYS'}

proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR',
                   'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']

ww = { 'A': 0.5,'R': 1.81,'N': 0.85,'D': 3.64,'C': -0.02,
       'Q': 0.77,'E': 3.63,'G': 1.15,'H': 2.33,'I': -1.12,
       'L': -1.25,'K': 2.80,'M': -0.67,'F': -1.71,'P': 0.14,
       'S': 0.46,'T': 0.25,'W': -2.09,'Y': -0.71,'V': -0.46 }

ww2 = { 'A': 0.17,'R': 0.81,'N': 0.42,'D': 1.23,'C': -0.24,
       'Q': 0.58,'E': 2.02,'G': 0.01,'H': 0.96,'I': -0.31,
       'L': -0.56,'K': 0.99,'M': -0.23,'F': -1.13,'P': 0.45,
       'S': 0.13,'T': 0.14,'W': -1.85,'Y': -0.94 ,'V': 0.07}


def featurizer(file, chains, database, processes):
    try:
        filename = file
        read = []
        with open(filename) as myfile:
            for line in myfile:
                if line.startswith('MODEL'):
                    break
                else:
                    read.append(line)
    except FileNotFoundError:
        print ('Not in upload directory, I\' ll try in PDB database')
        print ('Downloading PDB...')
        pdbl = PDBList()
        url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="
        pdbid = url+file.split("_")[0]
        open(os.path.join(os.path.dirname(sys.argv[0]), file + ".pdb"), "wb" ).write(urllib.request.urlopen(pdbid).read())
        filename = os.path.join(os.path.dirname(sys.argv[0]), file + '.pdb')
    
    try:
        mol = Molecule(filename)
    except RuntimeError:
        print ('Nah, neither in upload directory')
        sys.exit("Give a valid pdb, or upload a pdb file")  
        
    file = os.path.splitext(os.path.basename(file))[0]
    mol.dropFrames(keep = 0)
    if chains:
        print ('Selecting Chains...')
        chains_sel = ""
        for chain in chains[0]:
            chains_sel = chains_sel + chain + " "
        print (chains_sel)
        mol.filter('protein and chain ' + chains_sel)
        print ('protein and chain ' + chains_sel)
        file = file + '_chainA_fixed.pdb'
        filename2 = os.path.join(os.path.dirname(sys.argv[0]), 'outputs/prepared/fixed/' + file)
        filename = os.path.join(os.path.dirname(sys.argv[0]), 'outputs/prepared/' + file)
    else:
        mol.filter('protein')
        filename2 = os.path.join(os.path.dirname(sys.argv[0]), 'outputs/prepared/fixed/' + file + '_fixed.pdb')
    
    filename3 = filename
    print ('Preparing Protein...')
    for resnm in mol.resname:
         if resnm not in proteinResidues:
             mol.resname[mol.resname == resnm] = substitutions[resnm]
    mol = proteinPrepare(mol, pH=7.0, returnDetails=False)
    mol.write(filename2)
    mol.renumberResidues()
    mol.center(loc=(0, 0, 0), sel='all')
    mol.write(filename)
    
    #For charges
    command = "pdb2pqr30 --with-ph=7.0 --ff=PARSE --ffout=AMBER --neutraln --neutralc " + filename + " " + os.path.join(os.path.dirname(sys.argv[0]), file + '.pqr')
    os.system(command)
    u = mda.Universe(os.path.join(os.path.dirname(sys.argv[0]), file + '.pqr'))
    
    filename = filename2
    
    mol = Molecule(filename2)
    for resnm in mol.resname:
         if resnm not in proteinResidues:
             mol.resname[mol.resname == resnm] = substitutions[resnm]
    mol.write(filename2)
    
    #ProtDCal
    #It takes some time, so run it in parallel as a subprocess
    command = "python " + os.path.join(os.path.dirname(sys.argv[0]), "ML/train_ProtDCal.py " + file + ' ' + filename + " &")
    p = subprocess.Popen(command, shell=True,
    stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
    
    #PARSE PDB
    p2 = PDBParser(QUIET=True)
    structure = p2.get_structure("X", filename)
    model = structure[0]
    
    #DSSP
    print ('Calculating DSSP...')
    dssp = DSSP(model, filename, dssp='dssp')
    dssp = pd.DataFrame.from_dict(dssp.property_dict).T
    dssp = dssp.rename({0: 'DSSP index', 1: 'Amino acid', 2: 'Secondary structure', 3: 'ASA', 
    4: 'Phi', 5: 'Psi', 6: 'TCO', 7: 'KAPPA', 8: 'ALPHA', 9: 'NH-->O_1_relidx', 10: 'NH-->O_1_energy', 
    11: 'O-->NH_1_relidx', 12: 'O-->NH_1_energy', 13: 'NH-->O_2_relidx', 14: 'NH-->O_2_energy',
    15: 'O-->NH_2_relidx', 16: 'O-->NH_2_energy'}, axis='columns')
    dssp['Secondary structure'].replace('-', 'C', inplace = True)
    
    #FREESASA
    print ('Calculating FREESASA...')
    structure = freesasa.Structure(filename)
    result = freesasa.calc(structure, freesasa.Parameters({'n-slices' : 1000}))
    residueAreas = result.residueAreas()
    frsasa = pd.DataFrame(columns = ['residueType', 'residueNumber', 'total', 'mainChain', 'sideChain', 'polar', 'apolar',
       'hasRelativeAreas', 'relativeTotal', 'relativeMainChain',
       'relativeSideChain', 'relativePolar', 'relativeApolar']).set_index('residueType')
    for ch in residueAreas:
        for res in residueAreas[ch]:
            df = pd.DataFrame(residueAreas[ch][res].__dict__.items()).T
            df.columns = df.iloc[0]
            df = df[1:].set_index(['residueType'])
            df = df.astype(float, errors='ignore').round(2)
            frsasa = frsasa.append(df)
    frsasa.rename(columns={'total': 'all_atoms_abs', 'sideChain': 'side_chain_abs', 'mainChain': 'main_chain_abs', 'polar': 'all_polar_abs', 'apolar': 'non_polar_abs'}, inplace=True)
    frsasa = frsasa.set_index(dssp.index)

    #MSMS
    print ('Calculating MSMS...')
    rd = ResidueDepth(model)
    appender = []
    for k in rd.property_keys:
        x = rd.property_dict[k]
        residue = k[1]
        resnum = residue[1]
        resdepth = x[0]
        cadepth = x[1]
        appender.append((resnum, resdepth, cadepth))
    msms = pd.DataFrame(appender, columns=['resnum', 'res_depth', 'ca_depth'])
    resnum = msms['resnum']
    resnum.index = dssp.index
    
    #MERGE DATAFRAMES
    dssp = dssp.loc[:, ['Amino acid','Secondary structure', 'Phi', 'Psi', 'TCO', 'KAPPA', 'ALPHA']]
    frsasa = frsasa.loc[:, ['all_atoms_abs', 'side_chain_abs', 'main_chain_abs', 'all_polar_abs', 'non_polar_abs']]
    msms = msms.loc[:, ['res_depth', 'ca_depth']]
    msms.index = dssp.index
    
    #Sequence profiling
    print ('Calculating Sequence profiling...')
    hhb = 'hhblits'
    hhbdb = os.path.join(database, "uniclust30_2018_08")
    
    seq = str().join(dssp.loc[:,'Amino acid'])
    seq_prof, prof = getSequenceProfile(seq, hhb, hhbdb, processes, niter=3)
    
    seq_prof = seq_prof[1:]
    seq_prof.index = dssp.index
    
    seq_prof2 = []
    for index, row in seq_prof.iterrows():
        for column in seq_prof.columns:
            if row.seq == column:
                print (float(seq_prof.at[index, column]))
                seq_prof2.append(float(seq_prof.at[index, column]))
    seq_prof2 = pd.DataFrame(seq_prof2, columns=['Conservation score'])
    seq_prof2.index = dssp.index
    
    #Squared flunctuations
    print ('Calculating Squared flunctuations...')
    prot = parsePDB(filename)
    calphas = prot.select('calpha')
    gnm = GNM('Protein')
    gnm.buildKirchhoff(calphas)
    gnm.calcModes()
    sqf = calcSqFlucts(gnm)
    
    anm = ANM('prot ANM analysis')
    anm.buildHessian(calphas)
    anm.calcModes()
    sqf2 = calcSqFlucts(anm)
    
    sqf = pd.DataFrame(sqf)
    sqf2 = pd.DataFrame(sqf2)
    sqf.index = dssp.index
    sqf2.index = dssp.index
    sqf = sqf.rename({0: 'GNM'}, axis='columns')
    sqf2 = sqf2.rename({0: 'ANM'}, axis='columns')
    sqf['GNM'] = (sqf['GNM'] - sqf['GNM'].min()) / (sqf['GNM'].max() - sqf['GNM'].min())
    sqf2['ANM'] = (sqf2['ANM'] - sqf2['ANM'].min()) / (sqf2['ANM'].max() - sqf2['ANM'].min())
    
    #Create charges
    print ('Calculating charges...')
    charges = pd.DataFrame({'Charges': u.residues.charges})
    charges.index = dssp.index
    
    #Wimley-White hydrophobicity scales
    print ('Calculating Wimley-White hydrophobicity scales...')    
    hydroph = []
    hydroph3 = []
    i = 0
    for residue in seq:
        if (residue == 'H') and (int(charges.reset_index().iloc[i]['Charges']) == 0):
            hydroph.append(0.11)
            hydroph3.append(0.17)
        elif (residue == 'E') and (int(charges.reset_index().iloc[i]['Charges']) == 0):
            hydroph.append(0.11)
            hydroph3.append(-0.01)
        elif (residue == 'D') and (int(charges.reset_index().iloc[i]['Charges']) == 0):
            hydroph.append(0.43)
            hydroph3.append(-0.07)
        else:
            hydroph.append(ww[residue])
            hydroph3.append(ww2[residue])
        i = i + 1
    
    hydroph = pd.DataFrame({'Hydrophobicity': hydroph})
    hydroph.index = dssp.index
    hydroph3 = pd.DataFrame({'Hydrophobicity3': hydroph3})
    hydroph3.index = dssp.index
    
    hydroph2 = hydroph.copy()
    hydroph2.loc[hydroph2["Hydrophobicity"] > 0] = 0
    hydroph2.rename({"Hydrophobicity": "Hydrophobicity2"}, axis='columns', inplace=True)
    hydroph4 = hydroph3.copy()
    hydroph4.loc[hydroph4["Hydrophobicity3"] > 0] = 0
    hydroph4.rename({"Hydrophobicity3": "Hydrophobicity4"}, axis='columns', inplace=True)
    
    #Radius of gyration
    print ('Calculating Radius of gyration...')
    rog = []
    for resid in range(0, len(resnum)):
        rog.append(u.select_atoms('protein and resid ' + str(resid)).radius_of_gyration())
    
    rog = pd.DataFrame({'Radius of gyration': rog})
    rog.index = dssp.index
    
    features = [resnum, dssp, frsasa, msms, sqf, sqf2, seq_prof2, charges, hydroph, hydroph2, hydroph3, hydroph4, rog]
    features2 = pd.concat(features, axis=1)
    
    #Mean values based on nearby residues
    neigh_AA = ["neigh_A", "neigh_R", "neigh_N", "neigh_D", "neigh_C", "neigh_E", "neigh_Q", "neigh_G", "neigh_H", "neigh_I", "neigh_L", "neigh_K", "neigh_M",
          "neigh_F", "neigh_P", "neigh_S", "neigh_T", "neigh_W", "neigh_Y", "neigh_V"]
    
    print ('Calculating Mean values based on nearby residues...')
    cols = features2.drop(columns = ['Amino acid', 'Secondary structure', 'resnum']).columns
    ress = pd.DataFrame(resnum, resnum.index).reset_index()
    ress['resnum'] = list(range(0, len(resnum)))
    i = 0
    for ch in ress.level_0.unique():
        ress.loc[ress.level_0 == ch, 'level_0'] = i
        i = i + 1
    
    features2['resnum2'] = ress['resnum'].values
    prot = mda.Universe(filename3)
    residues2 = prot.select_atoms('protein and name CA').positions
    feature = pd.DataFrame(index = dssp.index)
    feature = feature.reindex(feature.columns.tolist() + neigh_AA, axis=1)
    global neighbors #nasty
    def neighbors(kkk):
        for index, resid in ress.iloc[kkk].iterrows():
            residues1 = prot.select_atoms('protein and resid ' + str(resid.resnum) + ' and name CA').positions
            dists3 = distance_array(residues1, residues2)
            dists3 = pd.DataFrame(dists3).T
            dists3.index = dssp.index
            distindex = dists3[dists3 < 7].dropna()
    
            neigh = features2.loc[distindex.index].groupby(features2.loc[distindex.index]['Amino acid'])['resnum'].nunique()
            for aa in ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]:
                if aa in neigh:
                    feature.loc[features2["resnum2"] == resid.resnum, 'neigh_' + aa] = neigh[aa]
                else:
                    feature.loc[features2["resnum2"] == resid.resnum, 'neigh_' + aa] = 0
    
            for column in cols:
                feature.loc[features2["resnum2"] == resid.resnum, 'mean_' + column] = features2.loc[distindex.index][column].mean()
            feature.loc[features2["resnum2"] == resid.resnum, '#Neighbors'] = len(features2.loc[distindex.index])
        return feature
    
    n_cpus = int(np.ceil((len(resnum)) / np.ceil((len(resnum)) / processes)))
    
    chunks = [range(0, len(resnum))[i::n_cpus] for i in range(n_cpus)]
    pool = Pool(processes=n_cpus)
    bla = pool.map(neighbors, chunks)
    pool.close()
    pool.join()
    
    features2 = pd.concat([features2, bla[0]], axis = 1)
    for bla2 in bla:
        features2 = features2.fillna(bla2)
    
    features2 = features2.drop(columns = ['resnum2'])
    
    features2.to_csv(os.path.join(os.path.dirname(sys.argv[0]), "outputs/features/" + file + ".csv"))
    
    p.communicate()
    p.wait()

    #Clean leftovers
# =============================================================================
#     command = "rm -v !('dreamm.py')"
#     os.system(command)
# =============================================================================
