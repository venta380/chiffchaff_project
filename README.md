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
# Population genetic analysis 
# chiffchaf_popgen.py
python package chiffchaf_popgen.py uses python packages like Pandas, numpy and Biopython to calculate some population genetic statistics. Please make sure you check the dependencies in the beginning of the file. 
# Dxy (absolute divergence between two populations)

maf(minor allele frequency) files can be made using angsd using below command. 

```
angsd -b $outdir/lists/Abeitinus -anc $ref -setMinDepth 5 -minInd 7 -nThreads 8 -out $outdir/Abeitinus -doMajorMinor 1 -doMaf 1 -gl 2 
angsd -b $outdir/lists/Tristis -anc $ref -setMinDepth 5 -minInd 7 -nThreads 8 -out $outdir/Tristis -doMajorMinor 1 -doMaf 1 -gl 2 
```

The function in  ANGSD_DXY_function in chiffchaf_popgen.py takes two "maf" files generated using ANGSD and a csv file with chromosome start, end and sites as input(example of format below).

```
,CHROM,BIN_START,BIN_END,Nsites,Nsites_SYM
0,Chr1,10000.0,20000.0,3542.0,3751.0
1,Chr1,20000.0,30000.0,6543.0,6578.0
2,Chr1,30000.0,40000.0,4837.0,4909.0
```

The ANGSD_DXY_function can be run using the below command. As the maf files are massive in their memory, please make sure you separate the maf files per each chromosome and run them separately. Mention the key word in the end 'ALLO' or 'SYMP' for allopatric or sympatric. I've used UPPMAX cluster with 16 cores with 512 GB ram. 


```
python dxy_angsd_all.py pop1.maf pop2.maf keys.csv output_name 'ALLO/SYMP' &
```
```
python dxy_angsd_all.py $file_path/Abeitinus_Chr1.csv           $file_path/Tristis_Chr1.csv             $file_path/keys_Chr1.csv Chr1_allo 'ALLO' &
python dxy_angsd_all.py $file_path/Abeitinus_Chr1A.csv          $file_path/Tristis_Chr1A.csv            $file_path/keys_Chr1A.csv Chr1A_allo 'ALLO' &
python dxy_angsd_all.py $file_path/Abeitinus_Chr2.csv           $file_path/Tristis_Chr2.csv             $file_path/keys_Chr2.csv Chr2_allo 'ALLO' &
python dxy_angsd_all.py $file_path/Abeitinus_Chr3.csv           $file_path/Tristis_Chr3.csv             $file_path/keys_Chr3.csv Chr3_allo 'ALLO' &
python dxy_angsd_all.py $file_path/Abeitinus_Chr4.csv           $file_path/Tristis_Chr4.csv             $file_path/keys_Chr4.csv Chr4_allo 'ALLO' &
python dxy_angsd_all.py $file_path/Abeitinus_Chr4A.csv          $file_path/Tristis_Chr4A.csv            $file_path/keys_Chr4A.csv Chr4A_allo 'ALLO' &
wait
```



