# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 12:04:59 2019

@author: alexis
"""

import sys, os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from joblib import load
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import fileinput

def predict(f):
    print ('Predicting protein-membrane interfaces')
    df = pd.read_csv(os.path.join(os.path.dirname(sys.argv[0]), '..', 'data/dataset_res_depth<2.5_selected.csv'), thousands = ',', low_memory=False).drop(columns=['Unnamed: 0'])
    df.apply(pd.to_numeric)
        
    unimport_feats_ProtDCal = pd.read_csv(os.path.join(os.path.dirname(sys.argv[0]), '..', 'data/unimport_feats.csv'), low_memory=False).drop(columns=['Unnamed: 0'])
    undf = pd.DataFrame(columns=sum(unimport_feats_ProtDCal.values.tolist(), []))
    
    intersection_cols = df.columns & undf.columns
    df = df.drop(intersection_cols, axis = 1) 
    
    df = pd.concat([df, undf], axis=1).fillna(0)
    df = df.reindex(sorted(df.columns), axis=1)
    
    standarscaler = StandardScaler()
    data_scaled_all = standarscaler.fit_transform(df.drop(columns=['class']))
    data_scaled_all = pd.DataFrame(data_scaled_all, columns = df.drop(columns=['class']).columns)
    
    #Read untested case
    dff2 = pd.read_csv(os.path.join(os.path.dirname(sys.argv[0]), 'outputs/features/' + f + '.csv'), thousands = ',', low_memory=False)
    tmpp = pd.read_csv(os.path.join(os.path.dirname(sys.argv[0]), 'ML/ProtDCal_v4.5/Outputs/' + f + '/' + f + '.csv'), thousands = ',', sep="\t", low_memory=False)
    tmpp = tmpp.drop_duplicates()
    if (len(dff2) != len(tmpp)) or (len(tmpp) == 0):
        sys.exit("Mismatch in ProtDCal")    
    
    tmpp['chain'] = tmpp['AA'].str.split('_').str[-3]
    tmpp['resnum'] = tmpp['AA'].str.split('_').str[-1].astype(int)
    tmpp = tmpp.sort_values(by=['chain', 'resnum']).reset_index().drop(columns=['index', 'chain', 'resnum'])
    dff2 = dff2.sort_values(by=['Unnamed: 0', 'resnum']).reset_index().drop(columns=['index'])
    
    
    dfff = pd.concat([dff2, tmpp], axis = 1)
    dfff = dfff.drop(columns=['AA'])
    
    chain_and_residues2 = dfff[['Unnamed: 0', 'resnum']]
    chain_and_residues2 = chain_and_residues2.reset_index().drop(columns=['index'])
    
    dfff = dfff[dfff['res_depth'] < 2.5] # Comment out to keep all residues
    dfff = dfff.reset_index()
    
    #Transform categorical features to numerical
    AA = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    AAdf = pd.DataFrame(columns=AA)
    for aa in AA:
        AAdf[aa] = np.where((aa == dfff['Amino acid']) , 1, np.nan)
    AAdf = AAdf.fillna(0)
    dfff = pd.concat([dfff, AAdf], axis = 1)
    
    SS = ["B", "C", "E", "G", "H", "I", "S", "T"]
    SSdf = pd.DataFrame(columns=SS)
    for ss in SS:
        SSdf[ss] = np.where((ss == dfff['Secondary structure']) , 1, np.nan)
    SSdf = SSdf.fillna(0)
    SS2 = ["Bss", "Css", "Ess", "Gss", "Hss", "Iss", "Sss", "Tss"]
    SSdf.rename(columns=dict(zip(SSdf.columns, SS2)), inplace=True)
    dfff = pd.concat([dfff, SSdf], axis = 1)
    aminoacids = dfff['Amino acid']
    dfff = dfff.drop(columns=['Amino acid', 'Secondary structure', 'Unnamed: 1', 'index'])
    
    dfff = dfff.replace({",":""}, regex=True)
    dfff.drop(columns=['Unnamed: 0']).apply(pd.to_numeric)
    dfff = dfff[dfff['res_depth'] < 2.5] # Comment out to keep all residues
    
    dfff2 = dfff
    
    loaded_model = load(os.path.join(os.path.dirname(sys.argv[0]), '..', 'data/final_classifier.sav'))
    
    
    chain_and_residues = dfff[['Unnamed: 0', 'resnum']]
    chain_and_residues= chain_and_residues.reset_index().drop(columns=['index'])
    
    intersection_cols = df.columns & dfff.columns
    dfff = dfff[intersection_cols]
    
    intersection_cols = dfff.columns & undf.columns
    dfff = dfff.drop(intersection_cols, axis = 1) 
    
    dfff = pd.concat([dfff, undf], axis=1).fillna(0)
    dfff = dfff.reindex(sorted(dfff.columns), axis=1)
    
    #Apply Scaling on the test set
    data_scaled_test = standarscaler.transform(dfff)
    data_scaled_test = pd.DataFrame(data_scaled_test, columns = dfff.columns)
    data_scaled_test = data_scaled_test[~data_scaled_test.isin([np.nan, np.inf, -np.inf]).any(1)]
    
    result = loaded_model.predict(data_scaled_test)
    result = pd.DataFrame(result, columns=['class'])
    results = pd.concat([chain_and_residues, aminoacids, result], axis = 1)
        
    if not results[results['class'] == 1]['resnum'].to_list():
        print ('Could not predict any residues inserting the membrane.')
    else:
        res = pd.DataFrame()
        AA = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "S", "K"]
        for aa in AA:
            res = pd.concat([results[(dfff2[aa] == 1) & (results['class'] == 1)], res])
            
        res.sort_values(inplace = True, by = ['Unnamed: 0', 'resnum'])
    
        res2 = results[results['class'] == 1][~results[results['class'] == 1].isin(res)].dropna()
        res2['resnum'] = res2['resnum'].astype(int)
        res2.sort_values(inplace = True, by = ['Unnamed: 0', 'resnum'])
    
        res = res[res['resnum'] >= 0]
        res2 = res2[res2['resnum'] >= 0]
    
        if not (res.empty or res2.empty):
            if f.endswith('_fixed.pdb'):
                filename = os.path.join(os.path.dirname(sys.argv[0]), 'outputs/prepared/fixed/' + f)
            else:
                filename = os.path.join(os.path.dirname(sys.argv[0]), 'outputs/prepared/fixed/' + f + '_fixed.pdb')
    
            with fileinput.FileInput(filename, inplace=True) as file:
                for line in file:
                    if (len(line) > 72) and (line[:72] != ''):
                        line = line[:72] + "    " + line[76:]
                        print(line.rstrip('\n'))
                    else:
                        print(line.rstrip('\n'))

            u = mda.Universe(filename)
            prot = u.select_atoms('protein')
    
            segid = 0
            for chain in results['Unnamed: 0'].unique():
                if segid > 0:
                    if ('residues' in locals()) and not (res['resnum'][res['Unnamed: 0'] == chain].empty):
                        residues = residues + prot.select_atoms('segid ' + np.unique(prot.segids)[segid] + ' and resid ' + 
                                res['resnum'][res['Unnamed: 0'] == chain].to_string(index=False, header = False).strip().replace('\n','')).residues
                    else:
                        if not res['resnum'][res['Unnamed: 0'] == chain].empty:
                            residues = prot.select_atoms('segid ' + np.unique(prot.segids)[segid] + ' and resid ' + 
                                    res['resnum'][res['Unnamed: 0'] == chain].to_string(index=False, header = False).strip().replace('\n','')).residues
                    if ('residues2' in locals()) and not (res2['resnum'][res2['Unnamed: 0'] == chain].empty):
                        residues2 = residues2 + prot.select_atoms('segid ' + np.unique(prot.segids)[segid] + ' and resid ' + 
                             res2['resnum'][res2['Unnamed: 0'] == chain].to_string(index=False, header = False).strip().replace('\n','')).residues
                    else:
                        if not res2['resnum'][res2['Unnamed: 0'] == chain].empty:
                            residues2 = prot.select_atoms('segid ' + np.unique(prot.segids)[segid] + ' and resid ' + 
                              res2['resnum'][res2['Unnamed: 0'] == chain].to_string(index=False, header = False).strip().replace('\n','')).residues
                   
                else:
                    if not res['resnum'][res['Unnamed: 0'] == chain].empty:
                        residues = prot.select_atoms('segid ' + np.unique(prot.segids)[segid] + ' and resid ' + 
                                res['resnum'][res['Unnamed: 0'] == chain].to_string(index=False, header = False).strip().replace('\n','')).residues
                    if not res2['resnum'][res2['Unnamed: 0'] == chain].empty:
                        residues2 = prot.select_atoms('segid ' + np.unique(prot.segids)[segid] + ' and resid ' + 
                              res2['resnum'][res2['Unnamed: 0'] == chain].to_string(index=False, header = False).strip().replace('\n','')).residues
                segid = segid + 1
            
            residues1_com = np.array([resid.atoms.center_of_mass() for resid in residues])
            residues2_com = np.array([resid.atoms.center_of_mass() for resid in residues2])
            dists = distance_array(residues1_com, residues2_com)
            dists2 = distance_array(residues1_com, residues1_com)
        
            kr = 0
            for dist in dists:
                for dis in dist:
                    if dis < 14:
                       res2 = pd.concat([res2, res.iloc[[kr]]])
                       break
                kr = kr + 1
            
            res2.sort_values(inplace = True, by = ['Unnamed: 0', 'resnum'])
            res = results[results['class'] == 1][~results[results['class'] == 1].isin(res2)].dropna()
            
        brokens = []
        for i in res2[res2['class'] == 1][['Unnamed: 0']]['Unnamed: 0'].unique():
            L = chain_and_residues2[chain_and_residues2['Unnamed: 0'] == i]['resnum'].tolist()
            start, end = L[0] - 3, L[-1] + 3
            missing = sorted(set(range(start, end + 1)).difference(L))
            for j in res2[(res2['class'] == 1) & (res2['Unnamed: 0'] == i)]['resnum'].tolist():
                fl = 0
                for k in range(1,4):
                    if (j - k in missing) or (j + k in missing):
                        brokens.append(1)
                        fl = 1
                        break
                if fl == 0:
                    brokens.append(0)
        res2['broken_chain'] = brokens
    
        print ('The residues: \n', res2[res2['class'] == 1][['Unnamed: 0', 'Amino acid', 'resnum', 'broken_chain']].to_string(index=False, header = False).strip(), '\nare predicted to insert the membrane.')
    
        res2.rename(columns={"Unnamed: 0": "chain"}, inplace = True)
        res2[['chain', 'resnum', 'Amino acid', 'broken_chain']].to_csv(os.path.join(os.path.dirname(sys.argv[0]), "outputs/prepared/fixed/" + f).split('.')[0].split('_chainA')[0] + ".csv", index=False)
        if res2.empty:
            print ('Could not predict any residues inserting the membrane.')
    
