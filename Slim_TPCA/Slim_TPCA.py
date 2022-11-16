# coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance
from sklearn.metrics import roc_curve, auc
import random
from scipy.optimize import curve_fit
import seaborn as sns
from scipy.stats import *
from scipy import stats
import copy
import re
from sklearn.preprocessing import StandardScaler
import pylab as py 
import warnings
warnings.filterwarnings('ignore')

def preproc(table, ref_col=1):
    if ref_col == False:
        return table
    else:
        table_clean = copy.copy(table)
        for col in table_clean.columns[ref_col:]: 
            soluble_fraction = table_clean.loc[:,col] / table.iloc[:, ref_col]
            table_clean.loc[:,col] = soluble_fraction
        return table_clean

def dist(table, ref_col=1, method='cityblock'):
    table_clean = preproc(table, ref_col)
    table_values = tuple(table.iloc[:,1:].values)
    dist_matrix = pd.DataFrame(distance.cdist(table_values, table_values, metric=method), index=table_clean.iloc[:,0], columns = table_clean.iloc[:,0])
    dist_matrix.index.name = ''
    dist_matrix.columns.name = ''
    return round(dist_matrix,6)

def pair_found(table, pair_table, ref_col=1):
    table_clean = preproc(table, ref_col)
    list_pro = list(table_clean.iloc[:,0])
    pair_table_found = pair_table[np.array([pair_table.iloc[:,0][i] in list_pro for i in range(len(pair_table))]) & 
                                  np.array([pair_table.iloc[:,1][i] in list_pro for i in range(len(pair_table))])].reset_index(drop=True)
    return pair_table_found
    
def roc(table, pair_table, ref_col=1, method='cityblock'):
    pair_table_found = pair_found(table, pair_table, ref_col)
    dist_matrix = dist(table, ref_col, method)
    roc_label, roc_score = [], []
    for i in range(len(pair_table_found)):
        pro_a, pro_b = pair_table_found.iloc[i,0], pair_table_found.iloc[i,1]
        roc_score.append(-1 * dist_matrix.loc[pro_a,pro_b]) 
        roc_label.append(1)
        dist_matrix.loc[pro_a, pro_b], dist_matrix.loc[pro_b, pro_a] = 0,0
    neg_values = np.triu(dist_matrix, k=0).flatten() 
    neg_values = neg_values[neg_values!=0]
    random.seed(42)
    roc_score += list(-1 * np.array(random.sample(list(neg_values),100000))) 
    roc_label += [0] * 100000
    fpr,tpr,threshold = roc_curve(roc_label, roc_score, pos_label = 1)
    roc_auc = auc(fpr,tpr)
    return fpr, tpr, round(roc_auc,4)

def roc_plot(table, pair_table, ref_col=1, method='cityblock'):
    fpr, tpr, roc_auc = roc(table, pair_table, ref_col, method)
    plt.figure(figsize=(3,3))
    plt.plot(fpr, tpr, label='AUC={}'.format(roc_auc))
    plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc=4, fontsize=12)
    plt.title('ROC Curve', fontsize=14)

def complex_found(table, complex_table, ref_col=1):
    table_clean = preproc(table, ref_col)
    complex_out = complex_table[complex_table['Organism']=='Human'].reset_index(drop=True)
    sub_found, num_sub_found = [], []
    for i in range(len(complex_out)):
        l_sub = complex_out['subunits(UniProt IDs)'][i].split(';')
        l_found = [sub for sub in l_sub if sub in list(table_clean.iloc[:,0])]
        sub_found.append(';'.join(l_found))
        num_sub_found.append(len(l_found))
    complex_out['Subunit_Found'] = sub_found
    complex_out['No_Subunit_Found'] = num_sub_found
    complex_out = complex_out[complex_out['No_Subunit_Found']>2].reset_index(drop=True)
    return complex_out

def complex_dist(table, complex_table, ref_col=1, method='cityblock'):
    complex_table_found = complex_found(table, complex_table, ref_col)
    dist_matrix = dist(table, ref_col, method)
    l_dist = []
    for i in range(len(complex_table_found)):
        l_sub = complex_table_found['Subunit_Found'][i].split(';')
        sub_matrix = dist_matrix.loc[l_sub, l_sub]
        avg_dist = np.nanmean(sub_matrix.replace(0, np.nan))
        l_dist.append(avg_dist)
    complex_table_found['Avg_Dist'] = l_dist
    complex_table_found['Avg_Dist_Derived'] = 1 / (1+np.array(l_dist))
    return complex_table_found

def random_n(table, complex_table, ref_col=1, method='cityblock', samplesize=10000):
    l_n = list(set(complex_found(table, complex_table, ref_col)['No_Subunit_Found'])) 
    l_n.sort()
    dic_out = {}
    pairs_dist_table = dist(table, ref_col, method)
    pairs_dist_table = pairs_dist_table.replace(0, np.nan)
    for num in l_n:
        l_n_dist = []
        random.seed(42)
        for i in range(samplesize):
            random_proteins = random.sample(list(pairs_dist_table.index), num)
            random_sub_table = pairs_dist_table.loc[random_proteins, random_proteins]
            l_n_dist.append(np.nanmean(random_sub_table))
        dic_out[num] = l_n_dist
    return dic_out

def complex_signature_sample(table, complex_table, ref_col=1, method='cityblock', samplesize=10000):
    complex_table_found = complex_dist(table, complex_table, ref_col)
    dic_random = random_n(table, complex_table, ref_col, method, samplesize)
    p_value, z_score = [], []
    for i in range(len(complex_table_found)):
        n = complex_table_found['No_Subunit_Found'][i]
        avg_dist = complex_table_found['Avg_Dist'][i]
        avg_dist_derived = complex_table_found['Avg_Dist_Derived'][i]
        l_random_n = dic_random[n]
        l_random_n_derived = list(1/(1+np.array(l_random_n))) + [avg_dist_derived]
        p_value.append(np.sum(np.array(l_random_n)<avg_dist)/samplesize)
        z_score.append(stats.zscore(l_random_n_derived)[-1])
    complex_table_found['TPCA_Sig_P-value'] = p_value
    complex_table_found['TPCA_Sig_Z-score'] = z_score
    random_table = pd.DataFrame(dic_random)
    return random_table, complex_table_found

def complex_signature_beta(table, complex_table, ref_col=1, method='cityblock', samplesize=500):
    complex_table_found = complex_dist(table, complex_table, ref_col)
    dic_random = random_n(table, complex_table, ref_col, method, samplesize)
    dic_beta = {}
    for k in dic_random:
        random_k = list(dic_random[k])
        random_k.sort()
        dic_beta[k] = getattr(stats, 'beta').fit(random_k[int(samplesize*0.025):-1*int(samplesize*0.025)]) 
    p_value = []
    for i in range(len(complex_table_found)):
        n = complex_table_found['No_Subunit_Found'][i]
        avg_dist = complex_table_found['Avg_Dist'][i]
        a, b, c, d = dic_beta[n]
        p_value.append(beta.cdf(avg_dist, a, b, c, d)) 
    complex_table_found['TPCA_Sig_P-value'] = p_value
    random_table = pd.DataFrame(dic_random)
    return random_table, complex_table_found

def align(table_1, table_2, ref_col=1):
    table_1_align = preproc(table_1, ref_col)
    table_2_align = preproc(table_2, ref_col)
    list_pro = list(set(table_1_align.iloc[:,0]) & set(table_2_align.iloc[:,0]))
    table_1_align = pd.merge(table_1_align, pd.DataFrame({table_1_align.columns[0]:list_pro}), on=table_1.columns[0])
    table_2_align = pd.merge(table_2_align, pd.DataFrame({table_2_align.columns[0]:list_pro}), on=table_2.columns[0])
    table_1_align = table_1_align.sort_values(by=table_1_align.columns[0]).reset_index(drop=True)
    table_2_align = table_2_align.sort_values(by=table_2_align.columns[0]).reset_index(drop=True)
    return table_1_align, table_2_align


def dynamic_complex_absolute_sample(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=10000):
    table_1_align, table_2_align = align(table_1, table_2, ref_col)
    complex_dist_1 = complex_dist(table_1_align, complex_table, ref_col, method)
    complex_dist_2 = complex_dist(table_2_align, complex_table, ref_col, method)
    complex_dist_change = pd.merge(complex_dist_1, complex_dist_2, on=list(complex_dist_1.columns[:-2]), suffixes=('_1','_2'))
    complex_dist_change['Avg_Dist_change'] = complex_dist_change['Avg_Dist_1'] - complex_dist_change['Avg_Dist_2']
    dic_random_1 = random_n(table_1_align, complex_table, ref_col, method, samplesize)
    dic_random_2 = random_n(table_2_align, complex_table, ref_col, method, samplesize)
    random_table = pd.DataFrame(dic_random_1) - pd.DataFrame(dic_random_2)
    p_value = []
    for i in range(len(complex_dist_change)):
        n = complex_dist_change['No_Subunit_Found'][i]
        avg_dist = complex_dist_change['Avg_Dist_change'][i]
        l_random_n = random_table[n]    
        p_value.append(np.sum(np.array(l_random_n)>avg_dist)/samplesize) 
    complex_dist_change['Dynamic_P'] = p_value
    return random_table, complex_dist_change

def dynamic_complex_relative_sample(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=10000):
    table_1_align, table_2_align = align(table_1, table_2, ref_col)
    complex_dist_1 = complex_dist(table_1_align, complex_table, ref_col, method)
    complex_dist_2 = complex_dist(table_2_align, complex_table, ref_col, method)
    complex_dist_change = pd.merge(complex_dist_1, complex_dist_2, on=list(complex_dist_1.columns[:-2]), suffixes=('_1','_2'))
    complex_dist_change['Avg_Dist_relative_change'] = (complex_dist_change['Avg_Dist_1'] - complex_dist_change['Avg_Dist_2']) / complex_dist_change['Avg_Dist_1']
    dic_random_1 = random_n(table_1_align, complex_table, ref_col, method, samplesize)
    dic_random_2 = random_n(table_2_align, complex_table, ref_col, method, samplesize)
    random_table = (pd.DataFrame(dic_random_1) - pd.DataFrame(dic_random_2)) / pd.DataFrame(dic_random_1) 
    p_value = []
    for i in range(len(complex_dist_change)):
        n = complex_dist_change['No_Subunit_Found'][i]
        avg_dist = complex_dist_change['Avg_Dist_relative_change'][i]
        l_random_n = random_table[n]
        p_value.append(np.sum(np.array(l_random_n)>avg_dist)/samplesize) 
    complex_dist_change['Dynamic_P'] = p_value    
    return random_table, complex_dist_change

def dynamic_complex_absolute_beta(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=500):
    table_1_align, table_2_align = align(table_1, table_2, ref_col)
    complex_dist_1 = complex_dist(table_1_align, complex_table, ref_col, method)
    complex_dist_2 = complex_dist(table_2_align, complex_table, ref_col, method)
    complex_dist_change = pd.merge(complex_dist_1, complex_dist_2, on=list(complex_dist_1.columns[:-2]), suffixes=('_1','_2'))
    complex_dist_change['Avg_Dist_change'] = complex_dist_change['Avg_Dist_1'] - complex_dist_change['Avg_Dist_2']
    dic_random_1 = random_n(table_1_align, complex_table, ref_col, method, samplesize)
    dic_random_2 = random_n(table_2_align, complex_table, ref_col, method, samplesize)
    random_table = pd.DataFrame(dic_random_1) - pd.DataFrame(dic_random_2)
    dic_beta = {}
    for k in random_table:
        random_k = list(random_table[k])
        random_k.sort()
        dic_beta[k] = getattr(stats, 'beta').fit(random_k[int(samplesize*0.025):-1*int(samplesize*0.025)])
    p_value = []
    for i in range(len(complex_dist_change)):
        n = complex_dist_change['No_Subunit_Found'][i]
        avg_dist = complex_dist_change['Avg_Dist_change'][i]
        a, b, c, d = dic_beta[n]
        p_value.append(beta.sf(avg_dist, a, b, c, d))
    complex_dist_change['Dynamic_P'] = p_value
    return complex_dist_change

def dynamic_complex_relative_beta(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=500):
    table_1_align, table_2_align = align(table_1, table_2, ref_col)
    complex_dist_1 = complex_dist(table_1_align, complex_table, ref_col, method)
    complex_dist_2 = complex_dist(table_2_align, complex_table, ref_col, method)
    complex_dist_change = pd.merge(complex_dist_1, complex_dist_2, on=list(complex_dist_1.columns[:-2]), suffixes=('_1','_2'))

    complex_dist_change['Avg_Dist_relative_change'] = (complex_dist_change['Avg_Dist_1'] - complex_dist_change['Avg_Dist_2']) / complex_dist_change['Avg_Dist_1']
    dic_random_1 = random_n(table_1_align, complex_table, ref_col, method, samplesize)
    dic_random_2 = random_n(table_2_align, complex_table, ref_col, method, samplesize)
    random_table = (pd.DataFrame(dic_random_1) - pd.DataFrame(dic_random_2)) / pd.DataFrame(dic_random_1) 
    dic_beta = {}
    for k in random_table:
        random_k = list(random_table[k])
        random_k.sort()
        dic_beta[k] = getattr(stats, 'beta').fit(random_k[int(samplesize*0.025):-1*int(samplesize*0.025)])
    p_value = []
    for i in range(len(complex_dist_change)):
        n = complex_dist_change['No_Subunit_Found'][i]
        avg_dist = complex_dist_change['Avg_Dist_relative_change'][i]
        a, b, c, d = dic_beta[n]
        p_value.append(beta.sf(avg_dist, a, b, c, d))
    complex_dist_change['Dynamic_P'] = p_value    
    return complex_dist_change
