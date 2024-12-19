# Population structure analysis pipeline
A pipeline using SNP variant information to conduct population structure analysis. 

![image](https://github.com/user-attachments/assets/574bd86f-1035-492c-add0-92ca89268790)

## Dependent software/script
- [PLINK1.9](https://www.cog-genomics.org/plink2/)
- [VCF2PCACluster](https://github.com/hewm2008/VCF2PCACluster)
- [ADMIXTURE](https://dalexander.github.io/admixture/download.html)
- [FastTree](https://morgannprice.github.io/fasttree/)
- [VCFtools](https://vcftools.github.io/)
- bgzip (part of [htslib](https://github.com/samtools/htslib))
- [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)

## What to input
Just vcf file containing your population variant information in compressed format (vcf.gz).
> Since ADMIXTURE only accpets vcf file with integer format chromosome ID (e.g. 01 for Chr01), please remember to modify your vcf file.

## What the pipeline does
- SNP filtered all missing site
- SNP vcf convert to plink format
- SNP pruned according to LD block
- PCA cluster
- ADMIXTURE analysis
- Phylogenic tree construction

## What to output
- PCA eigen principle component file used for PCA plot
- Tree file in nwk format
- ADMIXTURE Q file and P file
