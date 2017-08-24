# PCA (principal component analysis)
The R script PCA_chiffchaff.R takes input as a VCF file. This can be manually edited at line number 7 in the script in the variable "vcf.fn". 
# Population structure
The R script structure_chif_chaf_LEA.R uses the bioconductor package LEA which takes the .geno file generated by plink using the command 
```
plink --allow-extra-chr --chr-set 33 --recode --out plink  --vcf ./LD_venkat/Abeitinus.recode.vcf 
```
convert the ped to geno using ped2geno in software sNMF 
```
ped2geno plink.ped 
```
# Reference genome:
The Chiffchaff  reference used for this study was not denovo assembled but obtained by mapping the paired end reads to ficedula flycatcher genome. It can be obtained by bellow samtools commands using the bam files from ENA (European Nucleotide Archive) ERS1811978, ERS1811977, ERS1811976, ERS1811975. 
```
samtools mpileup  -u -f ficAlb2.fa --no-BAQ --count-orphans --min-BQ 5 --bam-list bamlist.txt | bcftools call -c -  > $TMPDIR/consensus_real.vcf
vcfutils.pl vcf2fq -d 5 -D 8000 $TMPDIR/consensus_real.vcf | gzip >  $TMPDIR/consensus_real.fq.gz
seqtk seq -a $TMPDIR/consensus_real.fq.gz > consensus_real.fa
```
A list ambiguous sites in consensus_real.fa are obtained by the script get_ambigious_sites.py by the below command
```
python get_ambigious_sites.py consensus_real.fa > ambiguous_sites.txt
```
Ambiguities were then assigned to random allele using the script fix_ambigious_sites.py
```
python fix_ambigious_sites.py consensus_real.fa ambiguous_sites.txt > consensus_fix.fa
```
You can also download the reference from https://www.dropbox.com/s/fdlh2hv5hogib6q/consensus_fix.fa.gz?dl=0
# Population genetic analysis 
# chiffchaf_popgen.py
python package chiffchaf_popgen.py uses python packages like Pandas, numpy and Biopython to calculate some population genetic statistics. Please make sure you check the dependencies in the beginning of the file. 
# Fst (genetic differentiation) between populations
Fst was calculated by below commands using ANGSD
```
angsd -b $outdir/lists2/Abeitinus -anc $ref -setMinDepth 5 -minInd 7 -nThreads 16 -out $outdir/Abeitinus_Sym -dosaf 1 -gl 2  -doMaf 8 -doMajorMinor 2 -doCounts 1 -fold 1 
angsd -b $outdir/lists2/Tristis -anc $ref -setMinDepth 5 -minInd 7 -nThreads 16 -out $outdir/Tristis_Sym -dosaf 1 -gl 2  -doMaf 8 -doMajorMinor 2 -doCounts 1 -fold 1
```

realSFS  $outdir/Abeitinus_Sym.saf.idx $outdir/Tristis_Sym.saf.idx -P 8 > $outdir/Abeitinus2.Tristis2.ml
realSFS fst index $outdir/Abeitinus_Sym.saf.idx $outdir/Tristis_Sym.saf.idx \
-sfs $outdir/Abeitinus2.Tristis2.ml -fstout $outdir/Abeitinus2.Tristis2


# Dxy (absolute divergence) between two populations

maf(minor allele frequency) files can be made using angsd using below command. 

```
angsd -b $outdir/lists/Abeitinus -anc $ref -setMinDepth 5 -minInd 7 -nThreads 8 -out $outdir/Abeitinus -doMajorMinor 1 -doMaf 1 -gl 2 
angsd -b $outdir/lists/Tristis -anc $ref -setMinDepth 5 -minInd 7 -nThreads 8 -out $outdir/Tristis -doMajorMinor 1 -doMaf 1 -gl 2 
```

The function in ANGSD_DXY_function in chiffchaf_popgen.py takes two "maf" files generated using ANGSD and a csv file with chromosome start, end and sites as input(example of format below).

```
,CHROM,BIN_START,BIN_END,,Nsites,Nsites_SYM
0,Chr1,10000.0,20000.0,3542.0,3751.0
1,Chr1,20000.0,30000.0,6543.0,6578.0
2,Chr1,30000.0,40000.0,4837.0,4909.0
```

The ANGSD_DXY_function can be run using the below command. As the maf files are massive in their memory, please make sure you separate the maf files per each chromosome and run them separately. Mention the key word in the end 'ALLO' or 'SYMP' for allopatric or sympatric. 


```
python dxy_angsd_all.py pop1.maf pop2.maf keys.csv output_name 'ALLO/SYMP' &
```
```
python dxy_angsd_all.py $file_path/Abeitinus_Chr1.csv           $file_path/Tristis_Chr1.csv             $file_path/keys_Chr1.csv Chr1_allo 'ALLO' 
python dxy_angsd_all.py $file_path/Abeitinus_Chr1A.csv          $file_path/Tristis_Chr1A.csv            $file_path/keys_Chr1A.csv Chr1A_allo 'ALLO' 
python dxy_angsd_all.py $file_path/Abeitinus_Chr2.csv           $file_path/Tristis_Chr2.csv             $file_path/keys_Chr2.csv Chr2_allo 'ALLO' 
python dxy_angsd_all.py $file_path/Abeitinus_Chr3.csv           $file_path/Tristis_Chr3.csv             $file_path/keys_Chr3.csv Chr3_allo 'ALLO' 
python dxy_angsd_all.py $file_path/Abeitinus_Chr4.csv           $file_path/Tristis_Chr4.csv             $file_path/keys_Chr4.csv Chr4_allo 'ALLO' 
python dxy_angsd_all.py $file_path/Abeitinus_Chr4A.csv          $file_path/Tristis_Chr4A.csv            $file_path/keys_Chr4A.csv Chr4A_allo 'ALLO' 
```
# LD (linkage disequilibrium analysis)
Perpare the files for the analysis using the filtering schema. 

```
vcftools --vcf Abeitinus.recode.vcf       --recode --out ./LD_venkat/Abeitinus        --max-missing-count 3
plink --allow-extra-chr --chr-set 33 --recode --out ./LD_venkat/Abeitinus        --vcf ./LD_venkat/Abeitinus.recode.vcf      
plink --allow-extra-chr --chr-set 33 --file ./LD_venkat/Abeitinus       --out ./LD_venkat/Abeitinus       --r2
#to check the averege r2 for windows of 10 kb accross genome
plink --allow-extra-chr --chr-set 33 --file ./LD_venkat/Abeitinus --ld-window 2 --ld-window-kb 10 --ld-window-r2 0 --maf 0.2 --out ./LD_venkat/Abeitinus_window --r2
```
use the outputs as the input for LD.py
```
LD.py Abeitinus_window.ld
```




