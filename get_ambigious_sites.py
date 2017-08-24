import sys
import string
import collections

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

fasta=str(sys.argv[1])

def fasta_dict(fasta):
        fasta = SeqIO.parse(fasta,"fasta")
        seq_dict={}
        for record in fasta:
                seq_dict[record.id]=record.seq
        return seq_dict

lengths_dict=fasta_dict(fasta)

for i in lengths_dict.iterkeys():
	for k in range(0, len(str(lengths_dict[i]))):
		j = str(lengths_dict[i])[k]
		if j == 'R' or j == 'Y' or j == 'S' or j == 'W' or j == 'K' or j == 'M' or j == 'B' or j == 'D' or j == 'H' or j == 'V' :
			print i + "\t" + str(k+1) + "\t" + str(k+1) 


