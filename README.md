# PCA (principal component analysis)
The R script PCA_chiffchaff.R takes input as a VCF file. This can be manually eddited at line number 7 in the script in the variable "vcf.fn". 
# Population structure
The R script structure_chif_chaf_LEA.R uses the bioconductor pakage LEA which takes the .geno file generated by plink using the comand 
```
plink --allow-extra-chr --chr-set 33 --recode --out plink  --vcf ./LD_venkat/Abeitinus.recode.vcf 
```
convert the ped to geno using ped2geno in softwere sNMF 
```
ped2geno plink.ped 
```
# Population genetic analysis 
# chiffchaf_popgen.py
python pakge chiffchaf_popgen.py uses python pakages like Pandas, numpy and Biopython to callulate some population genetic statestics
# Dxy (absolute divergence between two populations)
The function in ANGSD_DXY_function in chiffchaf_popgen.py takes two "maf" files generated using ANGSD and a csv file with chromosome start, end and sites as input
   
