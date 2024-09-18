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
        expand(f"pop_stru/ADMIXTURE/{vcf_basename}.K{{K}}.Q", K=K_values),
        expand(f"pop_stru/ADMIXTURE/{vcf_basename}.K{{K}}.P", K=K_values),
        expand(f"pop_stru/ADMIXTURE/log{{K}}.out", K=K_values),
        f"pop_stru/phylo/{vcf_basename}.nwk",
        f"pop_stru/PCA/{vcf_basename}.eigenvec",
        f"pop_stru/PCA/{vcf_basename}.eigenval"

rule VCF2plinkmap:
    input:
        vcf_file=config["vcf"]
    output:
        mapfile=f"plink/{vcf_basename}.map"
    log:
        "logs/vcf2plink.log"
    shell:
        """
        bcftools view \
        -H {input} \
        | cut -f 1 \
        | uniq \
        | awk '{{print $0"\t"$0}}' \
        > {output.mapfile}
        """

rule VCF2plinkped:
    input:
        vcf_file=config["vcf"],
        mapfile=f"plink/{vcf_basename}.map"
    output:
        pedfile=f"plink/{vcf_basename}.ped"
    log:
        "logs/vcf2plink.log"
    shell:
        """
        vcftools \
        --vcf {input.vcf_file} \
        --plink \
        --chrom-map {input.mapfile} \
        --out {output.pedfile}
        """
		
rule PLINKmakebed:
    input:
        f"plink/{vcf_basename}.ped",
        f"plink/{vcf_basename}.map"
    output:
        f"plink/{vcf_basename}.bed",
        f"plink/{vcf_basename}.bim",
        f"plink/{vcf_basename}.fam",
        f"plink/{vcf_basename}.nosex"
    log:
        "logs/plinkmakebed.log"
    params:
        f"plink/{vcf_basename}"
    shell:
        """
        plink \
        --file {params} \
        --make-bed \
        --double-id \
        --allow-extra-chr \
        --out {params} \
        &> {log}
        """

rule PLINKprune:
    input:
        f"plink/{vcf_basename}.bed",
        f"plink/{vcf_basename}.bim",
        f"plink/{vcf_basename}.fam"
    output:
        f"plink/{vcf_basename}.prune.in",
        f"plink/{vcf_basename}.prune.out"
    log:
        "logs/plinkprune.log"
    params:
        file=f"plink/{vcf_basename}",
        window_size=config["indep_pairwise"]["window_size"],
        step_size=config["indep_pairwise"]["step_size"],
        r2=config["indep_pairwise"]["r2"]
    shell:
        """
        plink \
        --allow-extra-chr \
        -bfile {params.file} \
        --indep-pairwise {params.window_size} {params.step_size} {params.r2} \
        --out {params.file} \
        &> {log}
        """

rule PLINKpruneExtract:
    input:
        f"plink/{vcf_basename}.bed",
        f"plink/{vcf_basename}.bim",
        f"plink/{vcf_basename}.fam",
        f"plink/{vcf_basename}.prune.in",
        f"plink/{vcf_basename}.prune.out"
    output:
        f"plink/{vcf_basename}.prune.bed",
        f"plink/{vcf_basename}.prune.bim",
        f"plink/{vcf_basename}.prune.fam",
        f"plink/{vcf_basename}.prune.nosex"
    log:
        "logs/plinkpruneextract.log"
    params:
        f"plink/{vcf_basename}",
        f"plink/{vcf_basename}.prune"
    shell:
        """
        plink \
        --allow-extra-chr \
        -bfile {params[0]} \
        --make-bed \
        --extract {input[3]} \
        --out {params[1]} \
        &> {log}
        """

rule PLINKprune2vcf:
    input:
        f"plink/{vcf_basename}.prune.bed",
        f"plink/{vcf_basename}.prune.bim",
        f"plink/{vcf_basename}.prune.fam",
        f"plink/{vcf_basename}.prune.nosex"
    output:
        f"{vcf_basename}.prune.vcf"
    log:
        "logs/plinkprune2vcf.log"
    params:
        f"plink/{vcf_basename}",
        f"{vcf_basename}.prune"
    threads: 16
    shell:
        """
        plink \
        --allow-extra-chr \
        --bfile {params[0]} \
        --export vcf \
        --out {params[1]} \
        &> {log}
        """

rule CompressPrunevcf:
    input:
        f"{vcf_basename}.prune.vcf"
    output:
        f"{vcf_basename}.prune.vcf.gz"
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
        f"{vcf_basename}.prune.vcf.gz"
    output:
        f"pop_stru/phylo/{vcf_basename}.prune.min4.fasta"
    log:
        "logs/vcf2phylip.log"
    shell:
        """
        python script/vcf2phylip.py \
        -i {input} \
        --output-folder phylo/ \
        -r \
        -p \
        -f \
        &> {log}
        """

rule PCA:
    input:
        vcf_file=config["vcf"]
    output:
        f"pop_stru/PCA/{vcf_basename}.eigenvec",
        f"pop_stru/PCA/{vcf_basename}.eigenval"
    threads: 32
    log:
        "logs/PCA.log"
    params:
        f"pop_stru/PCA/{vcf_basename}"
    shell:
        """
        VCF2PCACluster -InVCF {input} \
        -OutPut {params} \
        -MAF 0.05 \
        &> {log}
        """

rule ADMIXTURE:
    input:
        f"plink/{vcf_basename}.prune.bed"
    output:
        Q_file=f"pop_stru/ADMIXTURE/{vcf_basename}.K{{K}}.Q",
        P_file=f"pop_stru/ADMIXTURE/{vcf_basename}.K{{K}}.P",
        log_file=f"pop_stru/ADMIXTURE/log{{K}}.out"
    params:
        K=lambda wildcards: wildcards.K
    log:
        f"pop_stru/ADMIXTURE/log{{K}}.out"
    shell:
        """
        admixture --cv {input} {params.K} | tee {log}
        """

rule Phylogenictree:
    input:
        f"pop_stru/phylo/{vcf_basename}.prune.min4.fasta"
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
