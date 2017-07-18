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
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import gzip


#/home/venkat/.pyenv/versions/2.7.6/lib/python2.7/site-packages


def run_vcftools_population(population, population_1_sites, population_1_list):
    subprocess.call(["vcftools --vcf passed_vqsr_temp.vcf  --keep "+str(population_1_list)+ " --positions "+str(population_1_sites)+" --recode --out "+"./"+str(population)], shell=True)
    subprocess.call(["vcftools --vcf "+"./"+str(population)+".recode.vcf --site-pi --out "+"./"+str(population)+".PI.sites" ], shell=True)
    subprocess.call(["vcftools --vcf "+"./"+str(population)+".recode.vcf"+ "\t" + "--freq --out "+"./"+str(population)], shell=True)
 
def run_vcftools_generate_files(population_1_list, population_2_list, population_1, population_2, common_sites):
    subprocess.call(["vcftools --vcf passed_vqsr_temp.vcf --positions "+common_sites+" --recode --out $TMPDIR/output"], shell=True)
    subprocess.call(["vcftools --vcf $TMPDIR/output.recode.vcf  --weir-fst-pop "+str(population_1_list)+" --weir-fst-pop "+str(population_2_list)+" --out "+"./"+str(population_1+'_'+population_2)], shell=True)
    subprocess.call(["vcftools --vcf $TMPDIR/output.recode.vcf  --weir-fst-pop "+str(population_1_list)+" --weir-fst-pop "+str(population_2_list)+" --fst-window-size 10000 --out "+"./"+str(population_1+'_'+population_2)+"_fst_10000"], shell=True)



def  assess_coverage_filter_7X_induv(list_1, ind):
        file_1 =open(list_1[0])
        file_2 =open(list_1[1])
        file_3 =open(list_1[2])
        file_4 =open(list_1[3])
        file_5 =open(list_1[4])
        file_6 =open(list_1[5])
        file_7 =open(list_1[6])
        file_8 =open(list_1[7])
        file_9 =open(list_1[8])
        file_10=open(list_1[9])
        for line_1 in file_1:
                scaffold=line_1.strip().split()[0]
                position=line_1.strip().split()[1]
                cov_1=line_1.strip().split()[2]
                cov_2=file_2.readline().strip().split()[2]
                cov_3=file_3.readline().strip().split()[2]
                cov_4=file_4.readline().strip().split()[2]
                cov_5=file_5.readline().strip().split()[2]
                cov_6=file_6.readline().strip().split()[2]
                cov_7=file_7.readline().strip().split()[2]
                cov_8=file_8.readline().strip().split()[2]
                cov_9=file_9.readline().strip().split()[2]
                cov_10=file_10.readline().strip().split()[2]
                cover_list=[int(cov_1),int(cov_2),int(cov_3),int(cov_4),int(cov_5),int(cov_6),int(cov_7),int(cov_8),int(cov_9),int(cov_10)]
                count=0
                for l in cover_list:
                        if l >= 7:
                                count+=1
                if count >= ind:
                        DESSISION="P"
                else:
                        DESSISION="F"
                yield [scaffold, int(position), DESSISION]


def print_window_file(output, gen_gen, D_calls):
	outFileName = output
	outFile = open(outFileName, "w")
	outFile.write("Scaffold,Start,End,Nsites"+str(output)+"\n")
	with open (D_calls, 'r') as D_cal_vcftools:
	        first_line=D_cal_vcftools.readline()
	        D_value_dict={}
	        for n in D_cal_vcftools:
	                scaffold=n.strip().split()[0]
	                values=n.strip().split()[1:]
	                values.append(0)
	                if scaffold not in D_value_dict.keys():
	                        D_value_dict[scaffold]=[]
	                        D_value_dict[scaffold].append(values)
	                else:
	                        D_value_dict[scaffold].append(values)
	for i in gen_gen:
	        scaffold=i[0]
	        position=i[1]
	        DESSISION=i[2]
	        if DESSISION=="P":
	                if scaffold in D_value_dict.keys():
	                        #print scaffold
	                        for j in range(0,len(D_value_dict[scaffold])):
	                                start=int(D_value_dict[scaffold][j][0])
	                                end=int(D_value_dict[scaffold][j][1])
	                                if position >= start and position <= end:
	                                        D_value_dict[scaffold][j][2]+=1
	scaffold_track=[]
	for scaffold in D_value_dict.keys():
	        if scaffold not in scaffold_track:
	                scaffold_track.append(scaffold)
	                print scaffold
	        for window in D_value_dict[scaffold]:
	                if int(window[2]) > 0:
	                        #outFile.write("Scaffold,Start,End,Nsites")
	                        outFile.write(str(scaffold)+","+str(window[0])+","+str(window[1])+","+str(window[2])+"\n")
	                        sys.stdout.flush()
	

def fasta_dict(fasta):
	fasta = SeqIO.parse(fasta,"fasta")
	seq_dict = {}
	for record in fasta:
		seq_dict[record.id]=record.seq
	return seq_dict


def join_raw_data_base(lists):
	for i in range(1,len(lists)):
		j=i-1
		if i == 1:
			new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END'], how=outer)
		else:
			new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END'], how=outer)
		nextone=new_2
	return nextone


def join_data_base(lists):
	for i in range(1,len(lists)):
		j=i-1
		if i == 1:
			new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END'], how='outer', sort='BIN_START')
			new_2 = new_2.rename(columns={'N_VARIANTS_x': 'N_VARIANTS_'+str(j), 'WEIGHTED_FST_x': 'WEIGHTED_FST_'+str(j), 'MEAN_FST_x': 'MEAN_FST_'+str(j), 'N_VARIANTS_y': 'N_VARIANTS_'+str(i), 'WEIGHTED_FST_y': 'WEIGHTED_FST_'+str(i), 'MEAN_FST_y': 'MEAN_FST_'+str(i) })
		else:
			new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END'], how='outer', sort='BIN_START')
			new_2 = new_2.rename(columns={'N_VARIANTS_x': 'N_VARIANTS_'+str(j), 'WEIGHTED_FST_x': 'WEIGHTED_FST_'+str(j), 'MEAN_FST_x': 'MEAN_FST_'+str(j), 'N_VARIANTS_y': 'N_VARIANTS_'+str(i), 'WEIGHTED_FST_y': 'WEIGHTED_FST_'+str(i), 'MEAN_FST_y': 'MEAN_FST_'+str(i), 'N_VARIANTS': 'N_VARIANTS_'+str(i), 'WEIGHTED_FST': 'WEIGHTED_FST_'+str(i), 'MEAN_FST': 'MEAN_FST_'+str(i) })
		nextone=new_2
	return nextone

def callcualte_windowed_pi_from_site_pi(data_frame, site_pi_table,population):
	#final_window_pi={('irish_juvernica_final_PI'):'Filtered_sites_IJ',('kazak_juvernica_final_PI'):'Filtered_sites_KJ',('kaz_sin_final_PI'):'Filtered_sites_KS',('spanish_reali_final_PI'):'Filtered_sites_SR',('spanish_sinapis_final_PI'):'Filtered_sites_Sp_S',('swe_sin_final_PI'):'Filtered_sites_Sw_S'}
	win_sites=final_window_pi[population]
	for row_index, row in data_frame.iterrows():
		sites=row[win_sites]
		lists_pi=site_pi_table[(site_pi_table['CHROM'] == row.CHROM)]
		list_acc=range(row.BIN_START, row.BIN_END)
		sum_pi=np.nansum(lists_pi[lists_pi['POS'].isin(list_acc)].PI)
		lists_pi=np.nansum(lists_pi[lists_pi['POS'].isin(list_acc)].PI)/sites
		yield [row_index, lists_pi]


def extract_goods(file_name):
	good_sites=[]
	bad=[]
	with open (file_name) as f:
		first_line=f.readline()
		for line in f:
		    if (len(line.split())>=7):
		    	freq_list = [float(line.split()[-1][2:]), float(line.split()[-2][2:]), float(line.split()[-3][2:])]
		    	if 0.0 in freq_list:
		    		freq_list.sort()
		    		good_sites.append(line.strip())
		        if 0.0 not in freq_list:
		        	for n in freq_list:
		        		if not math.isnan(n):
		        			#print freq_list
		        			#print line
		        			bad.append(line.strip())
		    elif (len(line.split())==6):
		    	good_sites.append(line.strip())
	return good_sites

def convert_to_DF(file_name, header):
	#header=['CHROM', 'POS', 'N_ALLELES', 'N_CHR','al_1','al_2']
	fina_sites=[]
	df = pandas.DataFrame(fina_sites, columns=header)
	good_sites = extract_goods(file_name)
	for n in good_sites:
		nn=n.split()
		nn[1]=int(nn[1])
		#print nn[0]+ "/t" + nn[1]
		if (len(nn)>=7) or (int(nn[2])>2):
			if float(nn[-1][2:]) == 0.0:
				del nn[6]
				if len(nn) == 6:
					fina_sites.append(nn)
	
			elif float(nn[-2][2:]) == 0.0:
				del nn[5]
				if len(nn) == 6:
					fina_sites.append(nn)
			elif float(nn[-3][2:]) == 0.0:
				del nn[4]
				if len(nn) == 6:
					fina_sites.append(nn)
		else:
			df2= pandas.DataFrame({header[0]: [nn[0]],header[1]: [nn[1]],header[2]: [nn[2]],header[3]: [nn[3]],header[4]: [nn[4]],header[5]: [nn[5]] })
			df = df.append(df2, ignore_index=True)
			#print len(df)
	df3 = pandas.DataFrame(fina_sites, columns=header)
	lis_df=[df3, df]
	df = pandas.concat(lis_df, ignore_index=True)
	return df

		
		
def shared_private(p1,p2):
	fixed_diff=0
	private_pop1=0
	private_pop2=0
	shared=0
	fixed_same=0
	if ((p1 in [1.0, 0.0]) and (p2 not in [1.0, 0.0])):
		private_pop2+=1
	elif ((p1 in [1.0]) and (p2 in [0.0])):
		fixed_diff+=1
	elif ((p2 in [1.0]) and (p1 in [0.0])):
		fixed_diff+=1
	elif (p1 not in [1.0, 0.0]) and (p2 in [1.0, 0.0]):
		private_pop1+=1
	elif (p1 not in [1.0, 0.0]) and (p2 not in [1.0, 0.0]):
		shared+=1
	elif (p1 in [1.0, 0.0]) and (p2 in [1.0, 0.0]) and (p1==p2):
		fixed_same+=1
	return [fixed_diff, private_pop1, private_pop2, shared, fixed_same]



def ANGSD_freq_to_DF(file_1):
    file_a =pandas.read_table(file_1, compression='gzip')
    file_a =file_a.rename(columns={'chromo': 'CHROM', 'position': 'POS','anc': 'REF' ,'knownEM': 'minor_freq'}) 
    file_a =file_a[(file_a['REF']!='N')]
    file_a =file_a[['CHROM', 'POS', 'major','minor', 'minor_freq','nInd']]
    file_a['major_freq']= 1-file_a['minor_freq']
    return file_a



def fst_dxy(p1,p2):
	dxy=float((p1*(1-p2)) + (p2*(1-p1)))
	hs=((2.0*(p1*(1-p1)))+(2.0*(p2*(1-p2))))/2.0
	p=(p1+p2)/2.0
	q=((1-p1)+(1-p2))/2.0
	ht=2.0*p*q
	if ht>0.0:
		fst=((ht-hs))/ht
	else:
		fst=0
	return [fst,dxy]


def callculate_dxy_from_alle_frq_per_site(freq_file1,freq_file2):
	headera=['CHROM', 'POS', 'N_ALLELES', 'N_CHR','al_1_a','al_2_a']
	headerb=['CHROM', 'POS', 'N_ALLELES', 'N_CHR','al_1_b','al_2_b']
	a=convert_to_DF(freq_file1, headera)
	b=convert_to_DF(freq_file2, headerb)
	merged=pandas.merge(a,b, on=['CHROM','POS'])
	for row_index, row in merged.iterrows():
		fixed_pop1=0
		fixed_pop2=0
		private_pop1=0
		private_pop2=0
		shared=0
		if (row.al_1_a[0] == row.al_1_b[0]):
			p1=float(row.al_1_a[2:])
			p2=float(row.al_1_b[2:])
			dxy=float((p1*(1-p2)) + (p2*(1-p1)))
			pi_a=2*p1*(1-p1)
			pi_b=2*p2*(1-p2)
			if p1 in [1.0, 0.0] and p2 in [1.0, 0.0]:
				fixed_pop1+=1
				fixed_pop2+=1
			elif p1 in [1.0, 0.0] and p2 not in [1.0, 0.0]:
				fixed_pop1+=1
				private_pop2+=1
			elif p1 not in [1.0, 0.0] and p2 in [1.0, 0.0]:
				private_pop1+=1
				fixed_pop2+=1
			elif p1 not in [1.0, 0.0] and p2 not in [1.0, 0.0]:
				shared+=1
		elif (row.al_1_a[0] == row.al_2_b[0]):
			p1=float(row.al_1_a[2:])
			p2=float(row.al_1_b[2:])
			dxy=(p1*(1-p2)) + (p2*(1-p1))
			pi_a=2*p1*(1-p1)
			pi_b=2*p2*(1-p2)
			if p1 in [1.0, 0.0] and p2 in [1.0, 0.0]:
				fixed_pop1+=1
				fixed_pop2+=1
			elif p1 in [1.0, 0.0] and p2 not in [1.0, 0.0]:
				fixed_pop1+=1
				private_pop2+=1
			elif p1 not in [1.0, 0.0] and p2 in [1.0, 0.0]:
				private_pop1+=1
				fixed_pop2+=1
			elif p1 not in [1.0, 0.0] and p2 not in [1.0, 0.0]:
				shared+=1
		yield [row.CHROM, row.POS, dxy, pi_a, pi_b, fixed_pop1, private_pop1, fixed_pop2, private_pop2,shared]


#######################DXY_functions####################################

def dxy_window_function(file_a,file_b,nextone):
	headera=['CHROM', 'POS', 'N_ALLELES', 'N_CHR_a', 'al_1_', 'al_2_','al_3_','al_4_']
	headerb=['CHROM', 'POS', 'N_ALLELES', 'N_CHR_b', 'al_1_', 'al_2_','al_3_','al_4_']
	a=pandas.read_table(file_a, names=headera, engine='python',skiprows=1)
	b=pandas.read_table(file_b, names=headerb, engine='python',skiprows=1)
	merged=pandas.merge(a,b, on=['CHROM','POS'])
	merged=merged[(merged['N_CHR_a'] >= 14) &  (merged['N_CHR_b'] >= 14)]
	a=0
	b=0
	bads=[]
	goods=[]
	for window_index, window_1 in nextone.iterrows():
		window=merged[(merged['CHROM'] == window_1.CHROM)]
		window_fi=window[(window['POS']>=int(window_1.BIN_START)) & (window['POS']<=int(window_1.BIN_END))] 
		avg_dxy=[]
		avg_fst=[]
		fixed_diff=0
		private_pop1=0
		private_pop2=0
		shared=0
		fixed_same=0
		for row in window_fi.itertuples():
			allele_x={}
			allele_y={}
			if (int(row.N_ALLELES_x)<=2) or (int(row.N_ALLELES_y)<=2):
				#x alleles
				allele_x[row.al_1__x[0]]=float(row.al_1__x[2:])
				allele_x[row.al_2__x[0]]=float(row.al_2__x[2:])
				#y alleles
				allele_y[row.al_1__y[0]]=float(row.al_1__y[2:])
				allele_y[row.al_2__y[0]]=float(row.al_2__y[2:])
			elif (int(row.N_ALLELES_x)==3) or (int(row.N_ALLELES_y)==3):
				freq_list_x=[float(row.al_1__x[2:]), float(row.al_2__x[2:]),float(row.al_3__x[2:])]
				freq_list_y=[float(row.al_1__y[2:]), float(row.al_2__y[2:]),float(row.al_3__y[2:])]
				#x alleles
				if 0.0 in freq_list_x:
					if float(row.al_1__x[2:]) not in [0.0]: allele_x[row.al_1__x[0]]=float(row.al_1__x[2:])
					if float(row.al_2__x[2:]) not in [0.0]: allele_x[row.al_2__x[0]]=float(row.al_2__x[2:])
					if float(row.al_3__x[2:]) not in [0.0]: allele_x[row.al_3__x[0]]=float(row.al_3__x[2:])
				else:
					bads.append(row)
				#y alleles
				if 0.0 in freq_list_y:
					if float(row.al_1__y[2:]) not in [0.0]: allele_y[row.al_1__y[0]]=float(row.al_1__y[2:]) 
					if float(row.al_2__y[2:]) not in [0.0]: allele_y[row.al_2__y[0]]=float(row.al_2__y[2:])
					if float(row.al_3__y[2:]) not in [0.0]: allele_y[row.al_3__y[0]]=float(row.al_3__y[2:])
				else:
					bads.append(row)
			elif (int(row.N_ALLELES_x)==4) or (int(row.N_ALLELES_y)==4):
				freq_list_x=[float(row.al_1__x[2:]), float(row.al_2__x[2:]),float(row.al_3__x[2:])]
				freq_list_y=[float(row.al_1__y[2:]), float(row.al_2__y[2:]),float(row.al_3__y[2:])]
				#x alleles
				if 0.0 in freq_list_x:
					if float(row.al_1__x[2:]) not in [0.0]: allele_x[row.al_1__x[0]]=float(row.al_1__x[2:])
					if float(row.al_2__x[2:]) not in [0.0]: allele_x[row.al_2__x[0]]=float(row.al_2__x[2:])
					if float(row.al_3__x[2:]) not in [0.0]: allele_x[row.al_3__x[0]]=float(row.al_3__x[2:])
					if float(row.al_4__x[2:]) not in [0.0]: allele_x[row.al_4__x[0]]=float(row.al_4__x[2:])
				else:
					bads.append(row)
				#y alleles
				if 0.0 in freq_list_y:
					if float(row.al_1__y[2:]) not in [0.0]: allele_y[row.al_1__y[0]]=float(row.al_1__y[2:]) 
					if float(row.al_2__y[2:]) not in [0.0]: allele_y[row.al_2__y[0]]=float(row.al_2__y[2:])
					if float(row.al_3__y[2:]) not in [0.0]: allele_y[row.al_3__y[0]]=float(row.al_3__y[2:])
					if float(row.al_4__y[2:]) not in [0.0]: allele_y[row.al_4__y[0]]=float(row.al_4__y[2:])
				else:
					bads.append(row)
			if (len(allele_y) in [1,2]) and (len(allele_x) in [1,2]):
				x=list(allele_x.keys())
				y=list(allele_y.keys())
				if (set(x).issubset(set(y))) or (set(y).issubset(set(x))):
					if x[0] == y[0]:
						p1=allele_x[x[0]]
						p2=allele_y[y[0]]
					elif (len(x) ==1) and (x[0] == y[1]):
						p1=allele_x[x[0]]
						p2=allele_y[y[1]]
					elif (len(y) ==1) and (x[1] == y[0]):
						p1=allele_x[x[1]]
						p2=allele_y[y[0]]
					elif x[1] == y[0]:
						p1=allele_x[x[1]]
						p2=allele_y[y[0]]
					elif x[0] == y[1]:
						p1=allele_x[x[0]]
						p2=allele_y[y[1]]
					list_sp=shared_private(p1,p2)
					fixed_diff=fixed_diff+list_sp[0]
					private_pop1=private_pop1+list_sp[1]
					private_pop2=private_pop2+list_sp[2]
					shared=shared+list_sp[3]
					fixed_same=fixed_same+list_sp[4]
					dxy=float((p1*(1-p2)) + (p2*(1-p1)))
					hs=((2.0*(p1*(1-p1)))+(2.0*(p2*(1-p2))))/2.0
					p=(p1+p2)/2.0
					q=((1-p1)+(1-p2))/2.0
					ht=2.0*p*q
					if ht>0.0:
						fst=((ht-hs))/ht
					else:
						fst=0
					if dxy != 0.0:
						avg_fst.append(fst)
						avg_dxy.append(dxy)
		yield [window_index, np.nansum(avg_dxy)/(fixed_diff+private_pop1+private_pop2+shared), np.nanmean(avg_fst),fixed_diff, private_pop1, private_pop2, shared, fixed_same]




def ANGSD_DXY_function(file1,file2, fst_windows, sites):
    a=ANGSD_freq_to_DF(file1).round(3)
    b=ANGSD_freq_to_DF(file2).round(3)
    merge=pandas.merge(a,b, on=['CHROM','POS'])
    merge=merge[(merge['nInd_x'] >= 7) &  (merge['nInd_y'] >= 7)]
    for window_index, window_1 in fst_windows.iterrows():
		window=merge[(merge['CHROM'] == window_1.CHROM)]
		window_fi=window[(window['POS']>=int(window_1.BIN_START)) & (window['POS']<=int(window_1.BIN_END))]
		N_sites=window_1[sites]
		avg_dxy=[]
		avg_fst=[]
		avg_pi_a=[]
		avg_pi_b=[]
		fixed_diff=0
		private_pop1=0
		private_pop2=0
		shared=0
		fixed_same=0 
		for site_index, site in window_fi.iterrows():
			alleles_1={}
			alleles_1[site.major_x]=site.major_freq_x
			alleles_1[site.minor_x]=site.minor_freq_x
			alleles_2={}
			alleles_2[site.major_y]=site.major_freq_y
			alleles_2[site.minor_y]=site.minor_freq_y
			if site.major_x == site.major_y:
				p1 = site.major_freq_x
				p2 = site.major_freq_y
			elif site.major_x == site.minor_y:
				p1 = site.major_freq_x
				p2 = site.minor_freq_y
			elif site.minor_x == site.major_y:
				p1 = site.minor_freq_x
				p2 = site.major_freq_y
			dxy_TEMP=float((p1*(1-p2)) + (p2*(1-p1)))
			pi_a_temp=float(2*p1*(1-p1))
			pi_b_temp=float(2*p2*(1-p2))
			dxy=round(dxy_TEMP, 4)
			list_sp=shared_private(p1,p2)
			fixed_diff=fixed_diff+list_sp[0]
			private_pop1=private_pop1+list_sp[1]
			private_pop2=private_pop2+list_sp[2]
			shared=shared+list_sp[3]
			fixed_same=fixed_same+list_sp[4]
			hs=((2.0*(p1*(1-p1)))+(2.0*(p2*(1-p2))))/2.0
			p=(p1+p2)/2.0
			q=((1-p1)+(1-p2))/2.0
			ht=2.0*p*q
			if ht>0.0:
				fst=((ht-hs))/ht
			else:
				fst=0
			if dxy != 0.0:
				avg_fst.append(fst)
				avg_dxy.append(dxy)
				avg_pi_a.append(pi_a_temp)
				avg_pi_b.append(pi_b_temp)
		#if sorted(alleles_1.keys()) != sorted(alleles_2.keys()):
		#print [alleles_1.keys(),alleles_1.values(),alleles_2.keys(),alleles_2.values(), p1, p2, dxy, fst, pi_a, pi_b]
		yield [window_index, np.nansum(avg_dxy)/N_sites,np.nanmean(avg_fst),fixed_diff, private_pop1, private_pop2, shared, fixed_same, np.nansum(avg_pi_a)/N_sites, np.nansum(avg_pi_b)/N_sites]


def ANGSD_DXY_function_multiprocess(fst_windows, merge ,sites):
    for window_index, window_1 in fst_windows.iterrows():
		window=merge[(merge['CHROM'] == window_1.CHROM)]
		window_fi=window[(window['POS']>=int(window_1.BIN_START)) & (window['POS']<=int(window_1.BIN_END))]
		N_sites=window_1[sites]
		avg_dxy=[]
		avg_fst=[]
		avg_pi_a=[]
		avg_pi_b=[]
		fixed_diff=0
		private_pop1=0
		private_pop2=0
		shared=0
		fixed_same=0 
		for site_index, site in window_fi.iterrows():
			alleles_1={}
			alleles_1[site.major_x]=site.major_freq_x
			alleles_1[site.minor_x]=site.minor_freq_x
			alleles_2={}
			alleles_2[site.major_y]=site.major_freq_y
			alleles_2[site.minor_y]=site.minor_freq_y
			if site.major_x == site.major_y:
				p1 = site.major_freq_x
				p2 = site.major_freq_y
			elif site.major_x == site.minor_y:
				p1 = site.major_freq_x
				p2 = site.minor_freq_y
			elif site.minor_x == site.major_y:
				p1 = site.minor_freq_x
				p2 = site.major_freq_y
			dxy_TEMP=float((p1*(1-p2)) + (p2*(1-p1)))
			pi_a_temp=float(2*p1*(1-p1))
			pi_b_temp=float(2*p2*(1-p2))
			dxy=round(dxy_TEMP, 4)
			list_sp=shared_private(p1,p2)
			fixed_diff=fixed_diff+list_sp[0]
			private_pop1=private_pop1+list_sp[1]
			private_pop2=private_pop2+list_sp[2]
			shared=shared+list_sp[3]
			fixed_same=fixed_same+list_sp[4]
			hs=((2.0*(p1*(1-p1)))+(2.0*(p2*(1-p2))))/2.0
			p=(p1+p2)/2.0
			q=((1-p1)+(1-p2))/2.0
			ht=2.0*p*q
			if ht>0.0:
				fst=((ht-hs))/ht
			else:
				fst=0
			if dxy != 0.0:
				avg_fst.append(fst)
				avg_dxy.append(dxy)
				avg_pi_a.append(pi_a_temp)
				avg_pi_b.append(pi_b_temp)
		#if sorted(alleles_1.keys()) != sorted(alleles_2.keys()):
		#print [alleles_1.keys(),alleles_1.values(),alleles_2.keys(),alleles_2.values(), p1, p2, dxy, fst, pi_a, pi_b]
		print [window_index, np.nansum(avg_dxy)/N_sites,np.nanmean(avg_fst),fixed_diff, private_pop1, private_pop2, shared, fixed_same, np.nansum(avg_pi_a)/N_sites, np.nansum(avg_pi_b)/N_sites]



#######################END_of_DXY_functions####################################

