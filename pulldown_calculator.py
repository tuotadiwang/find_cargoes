#!/usr/bin/python

proteins = {'bait1':'pray1', 'bait2':'pray2'} # dictionary of bait_protein_names : pray_protein_names
bait_weights = {'bait1':80617.73, 'bait2':81802} # dictionary of bait_protein_names : molecular weight in Da
bait_concentrations = {'bait1':5.7, 'bait2':8.02} # dictionary of bait_protein_names : molecular concentration in mg/ml
pray_weights = {'pray1':15549.81, 'pray2':15549.81} # dictionary of pray_protein_names : molecular weight in Da
pray_concentrations = {'pray1':4.65, 'pray2':4.65} # dictionary of pray_protein_names : molecular concentration in mg/ml

mol_excess = 5
input_vol = 200 # uL
input_concentration = 0.0833 # unit: mg/ml
bait_amount = 50 # unit: ug
SDS = 40 # uL
DTT = 3 # uL
pulldown_vol = 500 # uL
slurry_vol = 60 # uL

import pandas as pd

def dict_drop_dup(dict):
    new_dict = {}
    for key,value in dict.items():
        if key not in new_dict.keys():
            new_dict[key] = value
    return new_dict

bait_mol = {}
pray_mol = {}
pray_amount = {}
bait_input = {}
pray_input = {}
bait_vol = {}
pray_vol = {}

for bait, pray in proteins.items():
    bait_mol[bait] = 50/bait_weights[bait]
    pray_mol[pray] = mol_excess * bait_mol[bait]
    pray_amount[pray] = pray_mol[pray] * pray_weights[pray]
    bait_input[bait] = input_vol * input_concentration / bait_concentrations[bait]
    pray_input[pray] = input_vol * input_concentration / pray_concentrations[pray]
    bait_vol[bait] = bait_amount / bait_concentrations[bait]
    pray_vol[pray] = pray_amount[pray] / pray_concentrations[pray]

inputs = {**dict_drop_dup(bait_input), **dict_drop_dup(pray_input)}
inputs_df = pd.DataFrame(inputs.items(), columns=['protein_name','protein'])
inputs_df['SDS_buffer'] = SDS
inputs_df['DTT'] = DTT
inputs_df['buffer'] = input_vol - SDS - DTT - inputs_df['protein']

pulldowns = pd.DataFrame(proteins.items(), columns=['bait_name','pray_name'])
pulldowns['bait'] = bait_vol.values()
pulldowns['pray'] = pray_vol.values()
pulldowns['slurry'] = slurry_vol
pulldowns['DTT'] = DTT
pulldowns['buffer'] = pulldown_vol - pulldowns['bait'] - pulldowns['pray'] - pulldowns['slurry'] - pulldowns['DTT']

print(inputs_df)
print(pulldowns)

