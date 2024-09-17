import os

configfile: "pop_structure_config.yaml"

vcf_basename = os.path.basename(config["vcf"])

if vcf_basename.endswith(".vcf.gz"):
    vcf_basename = vcf_basename[:-7]
elif vcf_basename.endswith(".vcf"):
    vcf_basename = vcf_basename[:-4]

K_values = config["K_values"]

rule all:
    input:
        expand(f"{admixture_dir}/{vcf_basename}.K{{K}}.Q", K=K_values),
        expand(f"{admixture_dir}/{vcf_basename}.K{{K}}.P", K=K_values),
        expand(f"{admixture_dir}/log{{K}}.out", K=K_values),


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
	--out {output[0]} \
        &> {log}
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
        --out {output[0]} \
        &> {log}
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
        --out {output[0]} \
        &> {log}
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
        --out {output[0]} \
        &> {log}
        """

rule PLINKprune2vcf:
    input:
        f"plink/{vcf_basename}.pruned",
        f"plink/{vcf_basename}.pruned.bed",
        f"plink/{vcf_basename}.pruned.bim",
        f"plink/{vcf_basename}.pruned.fam",
        f"plink/{vcf_basename}.pruned.nosex"
    output:
        f"{vcf_basename}.pruned",
        f"{vcf_basename}.pruned.vcf"
    log:
        "logs/plinkprune2vcf.log"
    threads: 16
    shell:
        """
        plink \
        --bfile {input[0]} \
        --export vcf \
        --out {output[0]} \
        &> {log}
        """

rule CompressPrunedvcf:
    input:
        f"{vcf_basename}.pruned.vcf"
    output:
        f"{vcf_basename}.pruned.vcf.gz"
    log:
        "logs/bgzip.log"
    threads: 16
    shell:
        """
        bgzip \
        -@ {threads} \
        {input} \
        &> {log}
        """

rule VCF2phylip:
    input:
        f"{vcf_basename}.pruned.vcf.gz"
    output:
        directory("pop_stru/phylo/"),
        f"pop_stru/phylo/{vcf_basename}.pruned.min4.fasta"
    log:
        "logs/vcf2phylip.log"
    shell:
        """
        python script/vcf2phylip.py \
        -i {input} \
        --output-folder {output[0]} \
        -r \
        -p \
        -f \
        &> {log}
        """

rule PCA:
    input:
        vcf_file=config["vcf"]
    output:
        "pop_stru/PCA/PCA",
        "pop_stru/PCA/PCA.eigenvec"
    threads: 32
    log:
        "logs/PCA.log"
    shell:
        """
        VCF2PCACluster -InVCF {input} \
        -OutPut {output[0]} \
        -MAF 0.05 \
        &> {log}
        """

rule ADMIXTURE:
    input:
        f"plink/{vcf_basename}.pruned.bed"
    output:
        Q_file=f"{admixture_dir}/{vcf_basename}.K{{K}}.Q",
        P_file=f"{admixture_dir}/{vcf_basename}.K{{K}}.P",
        log_file=f"{admixture_dir}/log{{K}}.out"
    params:
        K=lambda wildcards: wildcards.K  # 从 wildcards 中提取 K 值
    log:
        "{output.log_file}"
    shell:
        """
        admixture --cv {input} {params.K} | tee {log}
        """

rule Phylogenictree:
    input:
        f"pop_stru/phylo/{vcf_basename}.pruned.min4.fasta"
    output:
        f"pop_stru/phylo/{vcf_basename}.nwk"
    log:
        "logs/fasttree.log"
    shell:
        """
        FastTree \
        -gtr \
        -nt \
        {input} \
        > {output} \
        &> {log}
        """
