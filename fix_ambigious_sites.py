import random
import sys
import string
import collections

from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
from Bio.Seq import MutableSeq

fasta= str(sys.argv[1])
infile =  str(sys.argv[2])
#table = [line.strip().split() for line in open(infile)]



def fasta_dict(fasta):
        fasta = SeqIO.parse(fasta,"fasta")
        seq_dict={}
        for record in fasta:
                seq_dict[record.id]=MutableSeq(str(record.seq))
        return seq_dict

lengths_dict=fasta_dict(fasta)


with open(infile) as f:
	for line in f:
		chrom=line.split()[0]
		pos=int(line.split()[1])-1
		bases=[]
		info=line.split()[5:9]
		base1=int(info[0].split(':')[1])
		base2=int(info[1].split(':')[1])
		base3=int(info[2].split(':')[1])
		base4=int(info[3].split(':')[1])
		if base1 != 0:
			bases.append(info[0].split(':')[0])
		if base2 != 0:
			bases.append(info[1].split(':')[0])
		if base3 != 0:
			bases.append(info[2].split(':')[0])
		if base4 != 0:
			bases.append(info[3].split(':')[0])
		bases_list=[base1, base2, base3, base4] 
		maxbase=bases_list.index(max(bases_list))
		if maxbase==0:
			base_in_the_pos = "A"
		if maxbase==1:
			base_in_the_pos = "C"
		if maxbase==2:
			base_in_the_pos = "G"
		if maxbase==3:
			base_in_the_pos = "T"
		#print base_in_the_pos
		lengths_dict[chrom][pos] = base_in_the_pos


for i in lengths_dict:
        output=lengths_dict[i]
        for j in range(0, len(str(output))):
        	K=output[j]
        	if K == 'R':
        		output[j]=random.choice('AG')
        	if K == 'Y':
        		output[j]=random.choice('CT')
        	if K == 'S':
        		output[j]=random.choice('GC')
        	if K == 'W':
        		output[j]=random.choice('AT')
        	if K == 'K':
        		output[j]=random.choice('GT')
        	if K == 'M':
        		output[j]=random.choice('AC')
        	if K == 'B':
        		output[j]=random.choice('CGT')
        	if K == 'D':
        		output[j]=random.choice('AGT')
        	if K == 'H':
        		output[j]=random.choice('ACT')
        	if K == 'V':
        		output[j]=random.choice('ACG')
        print ">" + str(i)
        print output

