from IntegratedGradients import *
import pandas as pd
from tensorflow import keras
import numpy as np
import csv
import argparse
from tensorflow.python.framework.ops import disable_eager_execution
disable_eager_execution()

def one_hot(seq):
    one_hot_seq = np.zeros((4, len(seq)))
    for i in range(len(seq)):
            if seq[i] == "A" or seq[i] == 'a':
                one_hot_seq[0][i] = 1
            if seq[i] == "C" or seq[i] == 'c':
                one_hot_seq[1][i] = 1
            if seq[i] == "G" or seq[i] == 'g':
                one_hot_seq[2][i] = 1
            if seq[i] == "T" or seq[i] == 't':
                one_hot_seq[3][i] = 1
    return one_hot_seq

def make_prediction(seq, i, k, model_mf, model_sub):
    one_hot_seq = one_hot(seq[i:i+k])
    mid = one_hot(seq[i+int(k/2)])
    one_hot_seq = np.expand_dims(one_hot_seq, axis=-1)
    one_hot_seq = np.expand_dims(one_hot_seq, axis=0)

    mf_pred = np.exp(model_mf.predict(one_hot_seq)[0][0])

    sub_pred = model_sub.predict(one_hot_seq)
    sub_pred = np.multiply(1-mid[:,0],sub_pred[0])
    sub_pred = sub_pred/np.sum(sub_pred)
    return mf_pred, sub_pred

def explain(inputs,k,model_mf):
    ig = integrated_gradients(model_mf)

    outputs = []
    for seq_id, seq in enumerate(inputs,start=1):
        if len(seq) >= k:
            for i in range(len(seq) - k+1):
                one_hot_seq = one_hot(seq[i:i+k])
                one_hot_seq = np.expand_dims(one_hot_seq, axis=-1)

                exp = ig.explain(one_hot_seq)
                exp = np.squeeze(exp,axis=2)
                exp = np.transpose(exp)

                explanation = []
                for c in exp:
                    for num in c:
                        if num != 0:
                            explanation.append(num)

                outputs.append([seq_id,seq,explanation])
                
    df = pd.DataFrame(outputs,columns=['sequence id','sequence','explanation'])
    return df

def predict(inputs, k, model_mf, model_sub):
    outputs = []
    for seq_id, seq in enumerate(inputs,start=1):
        if len(seq) >= k:
            for i in range(len(seq) - k+1):
                mf_pred, sub_pred = make_prediction(seq, i, k, model_mf, model_sub)
                kmer_data = [seq_id, i+1, i+k, i+int(k/2)+1, seq[i+int(k/2)], mf_pred, sub_pred[0], sub_pred[1], sub_pred[2], sub_pred[3]]
                outputs.append(kmer_data)

    df = pd.DataFrame(outputs,columns=['sequence id', 'k-mer start', 'k-mer end', 'middle nucleotide position','middle nucleotide', 'mutation frequency', 'N>A', 'N>C', 'N>G', 'N>T'])
    return df
