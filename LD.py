import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import itertools
import math
import time
import sys
import chiffchaf_popgen
from pandas.tools.plotting import scatter_matrix

infile=str(sys.argv[1])



pwd='/proj/b2013146/nobackup/analysis/Faheema_analyis/vcftools/VCF_files/LD_venkat'
os.chdir(pwd)
data=pandas.read_csv('/proj/b2013146/nobackup/analysis/Faheema_analyis/vcftools/VCF_files/final_fst_dxy_FIXED_merged_db')[['CHROM', 'BIN_START','BIN_END']]


LD=pandas.read_table(infile, header=None,engine='python',skiprows=1,sep='\s+')
LD=LD.rename(columns={0: "CHROM_A",1: "START", 2: "SNP_A",3: "CHROM_B",4: "END", 5: "SNP_B",6: "R2"})


a=LD[(LD['CHROM_A'].isin(['ChrFAL35','ChrLGE22','Chr1A','Chr4A','ChrZ',"ChrFAL34"])==True)]
b=LD[(LD['CHROM_A'].isin(['ChrFAL35','ChrLGE22','Chr1A','Chr4A','ChrZ',"ChrFAL34"])==False)]
b['CHROM_A']="Chr"+b['CHROM_A']

frames = [a,b]
LD=pandas.concat(frames)



LD[['START', 'END','R2']]=LD[['START', 'END','R2']].apply(pandas.to_numeric)
data[['BIN_START', 'BIN_END',]]=data[['BIN_START', 'BIN_END']].apply(pandas.to_numeric)

LD=LD[(LD.CHROM_A == LD.CHROM_B)]
LD['DIST']=pandas.to_numeric(LD.END)-pandas.to_numeric(LD.START)
LD['R2_scaled']=LD['R2']/LD['DIST']


LD2=LD[(LD['DIST'] >= 10000) & (LD['DIST'] <= 20000)]
LD3=LD[(LD['DIST'] >= 20000) & (LD['DIST'] <= 30000)]



data['R2_scaled']=0.0
for window_index, window_1 in data.iterrows():
	window=LD[(LD.CHROM_A == window_1.CHROM) & (LD.START >= window_1.BIN_START) & (LD.END <= window_1.BIN_END)]
	if window.empty:
		data['R2_scaled'].set_value(window_index, value=0.0)
	else:
		data['R2_scaled'].set_value(window_index, value=np.nanmean(window.R2_scaled))


data.to_csv(infile+'final.ld')
