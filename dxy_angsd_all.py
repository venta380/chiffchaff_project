import sys
import os
import string
import pandas
import chiffchaf_popgen
import itertools
import subprocess
import sys
import string
from Bio import SeqIO
from Bio.Seq import Seq



pwd='/proj/b2013146/nobackup/analysis/Faheema_analyis/vcftools/VCF_files/temp_files/'
os.chdir(pwd)

#python dxy_angsd_all.py Abeitinus_file Tristis_file keys_file output_prefix 'ALLO'/'SYMP'


Abeitinus_maf=str(sys.argv[1])
Tristis_maf=str(sys.argv[2])
keys=str(sys.argv[3])
out_put=str(sys.argv[4])
style=str(sys.argv[5])


if style == 'ALLO':
	freq_list=['Abeitinus','Tristis']
elif style == 'SYMP':
	freq_list=['Abeitinus_Sym','Tristis_Sym']


primary_keys=pandas.read_csv(keys)


dxy_comb={str(i[0]+'_'+i[1]+'_dxy'): [pwd+(i[0]+'.frq'), (pwd+i[1]+'.frq')] for i in list(itertools.combinations(freq_list, 2))}
columns={str(i[0]+'_'+i[1]+'_dxy'): [(i[0]), (i[1])] for i in list(itertools.combinations(freq_list, 2))}
for i in dxy_comb.keys():
        filea=Abeitinus_maf
        fileb=Tristis_maf
        output_col=i
        fstcol=i[:-3]+'fst'
        fixedcol=i[:-3]+'fixed'
        privatea=i[:-3]+'private_a'
        privateb=i[:-3]+'private_b'
        sharedcol=i[:-3]+'shared'
        samecol=i[:-3]+'same_sites'
        pi_a=i[:-3]+'pi_a'
        pi_b=i[:-3]+'pi_b'
        primary_keys[str(output_col)]=0.0
        primary_keys[str(fstcol)]=0.0
        primary_keys[str(fixedcol)]=0.0
        primary_keys[str(privatea)]=0.0
        primary_keys[str(privateb)]=0.0
        primary_keys[str(sharedcol)]=0.0
        primary_keys[str(samecol)]=0.0
        primary_keys[str(pi_a)]=0.0
        primary_keys[str(pi_b)]=0.0
        for j in chiffchaf_popgen.ANGSD_DXY_function(filea,fileb,primary_keys, "Nsites"):
                #output_list [window_index, np.nanmean(avg_dxy), np.nanmean(avg_fst),fixed_diff, private_pop1, private_pop2, shared, fixed_same]
                primary_keys[str(output_col)]=primary_keys[str(output_col)].set_value(j[0], value=j[1])
                primary_keys[str(fstcol)]=primary_keys[str(fstcol)].set_value(j[0], value=j[2])
                primary_keys[str(fixedcol)]=primary_keys[str(fixedcol)].set_value(j[0], value=j[3])
                primary_keys[str(privatea)]=primary_keys[str(privatea)].set_value(j[0], value=j[4])
                primary_keys[str(privateb)]=primary_keys[str(privateb)].set_value(j[0], value=j[5])
                primary_keys[str(sharedcol)]=primary_keys[str(sharedcol)].set_value(j[0], value=j[6])
                primary_keys[str(samecol)]=primary_keys[str(samecol)].set_value(j[0], value=j[7])
                primary_keys[str(pi_a)]=primary_keys[str(pi_a)].set_value(j[0], value=j[8])
                primary_keys[str(pi_b)]=primary_keys[str(pi_b)].set_value(j[0], value=j[9])
                print j[0]
                sys.stdout.flush()

primary_keys.to_csv(pwd+out_put, )


