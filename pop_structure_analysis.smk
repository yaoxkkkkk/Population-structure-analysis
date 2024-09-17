import os

configfile: "pop_structure_config.yaml"

if vcf_basename.endswith(".vcf.gz"):
    vcf_basename = vcf_basename[:-7]
elif vcf_basename.endswith(".vcf"):
    vcf_basename = vcf_basename[:-4]

rule VCF2plink:
    input:
        vcf_file=config["vcf"]
    output:
        f"plink/{vcf_basename}",
        f"plink/{vcf_basename}.ped",
        f"plink/{vcf_basename}.map"
    log:
        "logs/vcf2plink.log"
    shell:
        """
        vcftools \
		    --vcf {input.vcf_file} \
		    --plink \
		    --out {output[0]}
        """

rule PLINKmakebed:
    input:
        f"plink/{vcf_basename}",
        f"plink/{vcf_basename}.ped",
        f"plink/{vcf_basename}.map"
    output:
        f"plink/{vcf_basename}",
        f"plink/{vcf_basename}.bed",
        f"plink/{vcf_basename}.bim",
        f"plink/{vcf_basename}.fam",
        f"plink/{vcf_basename}.nosex"
    log:
        "logs/plinkmakebed.log"
    shell:
        """
        plink \
        --file {input[0]} \
        --make-bed \
        --double-id \
        --out {output[0]}
        """

rule PLINKprune:
    input:
        f"plink/{vcf_basename}",
        f"plink/{vcf_basename}.bed",
        f"plink/{vcf_basename}.bim",
        f"plink/{vcf_basename}.fam"
    output:
        f"plink/{vcf_basename}.pruned",
        f"plink/{vcf_basename}.pruned.in",
        f"plink/{vcf_basename}.pruned.out"
    log:
        "logs/plinkprune.log"
    shell:
        """
        plink \
        -bfile {input[0]} \
        --indep-pairwise 500 50 0.2 \
        --out {output[0]}
        """

rule PLINKpruneExtract:
    input:
        f"plink/{vcf_basename}",
        f"plink/{vcf_basename}.pruned",
        f"plink/{vcf_basename}.bed",
        f"plink/{vcf_basename}.bim",
        f"plink/{vcf_basename}.fam",
        f"plink/{vcf_basename}.pruned.in",
        f"plink/{vcf_basename}.pruned.out"
    output:
        f"plink/{vcf_basename}.pruned",
        f"plink/{vcf_basename}.pruned.bed",
        f"plink/{vcf_basename}.pruned.bim",
        f"plink/{vcf_basename}.pruned.fam",
        f"plink/{vcf_basename}.pruned.nosex"
    log:
        "logs/plinkpruneextract.log"
    shell:
        """
        plink \
        -bfile {input[0]} \
        --make-bed \
        --extract {input[5]} \
        --out {output[0]}
        """

rule PCA_analysis:
    input:
        vcf_file=config["vcf"]
    output:
        "pop_stru/PCA/PCA",
        "pop_stru/PCA/PCA.eigenvec"
    threads: 32
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        VCF2PCACluster -InVCF {input} \
        -OutPut {output[0]} \
        -MAF 0.05
        """
rule 
