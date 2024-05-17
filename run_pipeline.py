import os
from Bio import SeqIO
import predict
from tensorflow import keras
import pandas as pd
import numpy as np
import h5py

# Assign 15-mer model
model_mf = keras.models.load_model('models/model_15_mf.h5')
model_sub = keras.models.load_model('models/model_15_sub.h5')

# Assign 15-mer dataset
f = h5py.File('data/IGHV_kmer15_by_family_cdr3frameshift.h5')

# IMGT gapped sequence for subregion assignment
seqs = SeqIO.to_dict(SeqIO.parse('data/imgt_vregion_gapped.fasta','fasta'))

# Motif of interest - AGCT 
OHS_motifs = ['AGCA','AGCC','AGCT','TGCA','TGCC','TGCT']
motif = 'AGCT'
dfs = []

# Across 7 families
for family in ['IGHV'+str(i) for i in range(1,8)]:

# Get available allele data
    alleles = f[family]['annotations']['allele'].tolist()
    alleles = [i.decode('UTF-8') for i in alleles]

    for allele in set(alleles):
        allele_indexes = [i for i in range(len(alleles)) if alleles[i]==allele]

        mid_bases = (f[family]['middle_base'][allele_indexes])
        mid_bases = [i.decode('UTF-8') for i in mid_bases]

        ind_bases = [(mid_bases[i],allele_indexes[i]) for i in range(len(mid_bases))] 

# Use only k-mers corresponding to motif identity
        motif_indexes = []
        for i in ind_bases:
            if i[0] == 'G':
                fivemer = (f[family]['fivemer'][i[1]].decode('UTF-8'))
                if fivemer[1]==motif[0] and fivemer[3]==motif[2] and fivemer[4]==motif[3]:
                    motif_indexes.append(i[1])
            if i[0] == 'C':
                fivemer = (f[family]['fivemer'][i[1]].decode('UTF-8'))
                if fivemer[0]==motif[0] and fivemer[1]==motif[1] and fivemer[3]==motif[3]:
                    motif_indexes.append(i[1])

# Get sequence and observed mutation frequency
        kmers = [f[family]['dna_kmer'][i].decode('UTF-8') for i in motif_indexes]
        freqs = [f[family]['freq'][i] for i in motif_indexes]

# Assign IMGT position and subregion
        positions = [f[family]['annotations']['offset'][i]+7 for i in motif_indexes]
        positions = [pos+seqs[[i for i in seqs.keys() if allele in i][0]].seq[:pos].count('.') for pos in positions]

        subregions = []
        for i in positions:
            if i <= 78: 
                subregions.append('FW1')
            elif i>78 and i<=114:
                subregions.append('CDR1')
            elif i>114 and i<=165:
                subregions.append('FW2')
            elif i>165 and i<=195:
                subregions.append('CDR2')
            elif i>195 and i<=312:
                subregions.append('FW3')
            elif i>312:
                subregions.append('CDR3')

# DeepSHM prediction
        preds = predict.predict(kmers,15,model_mf,model_sub)
        preds = preds['mutation frequency'].tolist()

# Integrated gradients explanation
        exp = predict.explain(kmers,15,model_mf)

# Format dataframe
        exp['allele'] = allele
        exp = exp.set_index('allele')

# Add metadata
        exp['Observed mutation frequency']=freqs
        exp['DeepSHM predicted mutation frequency']=preds
        exp['IMGT Position']=positions
        exp['IMGT Subregion']=subregions

        dfs.append(exp)

# Save as csv
df = pd.concat(dfs)
df.to_csv('data/'+motif+'_pos_truepred_igexps.csv')
