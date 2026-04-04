########################################
# CONFIG
########################################
configfile: "config.yaml"

SAMPLES = [s.strip() for s in open(config["sample_file"])]

OUT = config["outdir"]
FQ  = config["fq_dir"]

########################################
# PATHS
########################################
def fmt(p):
    return p.format(
        fq_dir=config["fq_dir"],
        outdir=config["outdir"],
        sample="{sample}"   # ✅ KEEP WILDCARD
    )

R1        = fmt(config["paths"]["R1"])
R2        = fmt(config["paths"]["R2"])

TRIM_R1   = fmt(config["paths"]["TRIM_R1"])
TRIM_R2   = fmt(config["paths"]["TRIM_R2"])

BAM_SORT  = fmt(config["paths"]["BAM_SORT"])
BAM_MD    = fmt(config["paths"]["BAM_MD"])
BAI       = fmt(config["paths"]["BAI"])

###########################################################
# -------- QC BLOCK and Trimming Data --------------------
###########################################################

# Raw Data FastQC
rule fastqc_raw:
    input: R1, R2
    output:
        f"{OUT}/pretrim_fastQC/{{sample}}_R1_fastqc.html",
        f"{OUT}/pretrim_fastQC/{{sample}}_R2_fastqc.html"
    threads: config["threads"]["fastqc"]
    shell:
        "mkdir -p {OUT}/pretrim_fastQC; fastqc -t {threads} -o {OUT}/pretrim_fastQC {input}"

# Raw Data MultiQC
rule multiqc_pre:
    input:
        expand(f"{OUT}/pretrim_fastQC/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        expand(f"{OUT}/pretrim_fastQC/{{sample}}_R2_fastqc.html", sample=SAMPLES)
    output:
        f"{OUT}/pretrim_fastQC/multiqc_report.html"
    shell:
        "multiqc {OUT}/pretrim_fastQC -o {OUT}/pretrim_fastQC -n multiqc_report.html"

# Trimming Data: TrimGalore
rule trim_galore:
    input: R1, R2
    output: TRIM_R1, TRIM_R2
    threads: config["threads"]["trim"]
    shell:
        "mkdir -p {OUT}/TrimGalore; trim_galore --paired --cores {threads} -o {OUT}/TrimGalore {input}"

# Trimmed Data FastQC
rule fastqc_trim:
    input: TRIM_R1, TRIM_R2
    output:
        f"{OUT}/post_trim_fastQC/{{sample}}_R1_val_1_fastqc.html",
        f"{OUT}/post_trim_fastQC/{{sample}}_R2_val_2_fastqc.html"
    threads: config["threads"]["fastqc"]
    shell:
        "mkdir -p {OUT}/post_trim_fastQC; fastqc -t {threads} -o {OUT}/post_trim_fastQC {input}"

# Trimmed Data MultiQC
rule multiqc_post:
    input:
        expand(f"{OUT}/post_trim_fastQC/{{sample}}_R1_val_1_fastqc.html", sample=SAMPLES),
        expand(f"{OUT}/post_trim_fastQC/{{sample}}_R2_val_2_fastqc.html", sample=SAMPLES)
    output:
        f"{OUT}/post_trim_fastQC/multiqc_report.html"
    shell:
        "multiqc {OUT}/post_trim_fastQC -o {OUT}/post_trim_fastQC -n multiqc_report.html"

########################################
# -------- ALIGNMENT BLOCK -------------
########################################

# Alignment: DragMap
rule align:
    input: TRIM_R1, TRIM_R2
    output: BAM_SORT
    threads: config["threads"]["dragmap"]
    params:
        ref=config["ref_genome"],
        samthreads=config["threads"]["samtools"]
    shell:
        """
        mkdir -p {OUT}/bam
        dragen-os -r {params.ref} \
            -1 {input[0]} -2 {input[1]} \
            --num-threads {threads} \
            --RGID illumina \
            --RGSM {wildcards.sample} | \
        samtools sort -@ {params.samthreads} -o {output}
        """

# MarkDup
rule markdup:
    input: BAM_SORT
    output:
        bam=BAM_MD,
        metrics=f"{OUT}/picard_info/{{sample}}_MarkDup_metrics.txt"
    threads: config["threads"]["picard"]
    shell:
        """
        mkdir -p {OUT}/picard_info
        picard MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            ASSUME_SORT_ORDER=coordinate
        """

# BAM Index
rule index:
    input: BAM_MD
    output: BAI
    threads: config["threads"]["samtools"]
    shell:
        "samtools index -@ {threads} {input} {output}"

# BAM QC
rule bamqc:
    input: bam=BAM_MD, bai=BAI
    output:
        file=f"{OUT}/bam/flagstat/{{sample}}.txt"
    threads: config["threads"]["samtools"]
    shell:
        """
        mkdir -p {OUT}/bam/flagstat
        samtools flagstat -@ {threads} {input.bam} > {output.file}
        """
# multiQC
rule multiqc_bam:
    input:
        expand(f"{OUT}/bam/flagstat/{{sample}}.txt", sample=SAMPLES)
    output:
        f"{OUT}/bam/flagstat/multiqc_report.html"
    shell:
        """
        multiqc {OUT}/bam/flagstat -o {OUT}/bam/flagstat -n multiqc_report.html
        """
        
########################################
# -------- SOMATIC ---------------------
########################################

# Mutect2 (tumor-only) -> generates unfiltered VCF + tumor.f1r2.tar.gz
rule mutect2:
    input: bam=BAM_MD, bai=BAI
    output:
        vcf=f"{OUT}/SomVar_VCF/{{sample}}.vcf.gz",
        stats=f"{OUT}/SomVar_VCF/{{sample}}.vcf.gz.stats",
        f1r2=f"{OUT}/SomVar_VCF/{{sample}}.f1r2.tar.gz"
    threads: config["threads"]["gatk"]
    params:
        ref=config["ref_fasta"],
        germ=config["germline_resource"],
        pon=config["pon"]
    shell:
        """
        mkdir -p {OUT}/SomVar_VCF
        gatk Mutect2 \
            -R {params.ref} \
            -I {input.bam} \
            --germline-resource {params.germ} \
            --panel-of-normals {params.pon} \
            --f1r2-tar-gz {output.f1r2} \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf}
        """

# Learn read-orientation model
rule learn_orientation:
    input:
        f"{OUT}/SomVar_VCF/{{sample}}.f1r2.tar.gz"
    output:
        f"{OUT}/SomVar_VCF/{{sample}}.artifact-priors.tar.gz"
    shell:
        "gatk LearnReadOrientationModel -I {input} -O {output}"

# Filter Mutect2 calls (apply filters + orientation priors)
rule filter:
    input:
        vcf=f"{OUT}/SomVar_VCF/{{sample}}.vcf.gz",
        stats=f"{OUT}/SomVar_VCF/{{sample}}.vcf.gz.stats",
        priors=f"{OUT}/SomVar_VCF/{{sample}}.artifact-priors.tar.gz"
    output:
        f"{OUT}/SomVar_VCF/{{sample}}.filtered.vcf.gz"
    params:
        ref=config["ref_fasta"]
    shell:
        """
        gatk FilterMutectCalls \
            -V {input.vcf} \
            --stats {input.stats} \
            --orientation-bias-artifact-priors {input.priors} \
            -R {params.ref} \
            -O {output}
        """

########################################
# -------- GERMLINE --------------------
########################################

# Build the DRAGEN STR Model from the Aligned Reads
rule dragstr_model:
    input: bam=BAM_MD, bai=BAI
    output:
        f"{OUT}/dragstr/{{sample}}_dragstr_model.txt"
    params:
        ref=config["ref_genome"]
    shell:
        """
        mkdir -p {OUT}/dragstr
        gatk CalibrateDragstrModel \
        -R {params.ref}/hg38.fa \
        -I {input.bam} \
        -str {params.ref}/str_table.tsv \
        -O {output}
        """

# HaplotypeCaller in DRAGEN mode
rule haplotypecaller:
    input:
        bam=BAM_MD,
        bai=BAI,
        drag=f"{OUT}/dragstr/{{sample}}_dragstr_model.txt"
    output:
        f"{OUT}/gvcf/{{sample}}.g.vcf.gz"
    threads: config["threads"]["gatk"]
    params:
        ref=config["ref_genome"]
    shell:
        """
        mkdir -p {OUT}/gvcf
        gatk --java-options '-Xmx12g' HaplotypeCaller \
            -R {params.ref}/hg38.fa \
            -I {input.bam} \
            -O {output} \
            -ERC GVCF \
            --dragen-mode true \
            --dragstr-params-path {input.drag} \
            --native-pair-hmm-threads {threads}
        """

########################################
# -------- TARGET RULES ----------------
########################################

rule qc:
    input:
        f"{OUT}/pretrim_fastQC/multiqc_report.html",
        f"{OUT}/post_trim_fastQC/multiqc_report.html"

rule bamqc_all:
    input:
        expand(f"{OUT}/bam/flagstat/{{sample}}.txt", sample=SAMPLES)

rule somatic:
    input:
        expand(f"{OUT}/SomVar_VCF/{{sample}}.filtered.vcf.gz", sample=SAMPLES)

rule germline:
    input:
        expand(f"{OUT}/gvcf/{{sample}}.g.vcf.gz", sample=SAMPLES)

rule somatic_full:
    input:
        expand(f"{OUT}/SomVar_VCF/{{sample}}.filtered.vcf.gz", sample=SAMPLES)

rule germline_full:
    input:
        expand(f"{OUT}/gvcf/{{sample}}.g.vcf.gz", sample=SAMPLES)

rule all:
    input:
        expand(f"{OUT}/gvcf/{{sample}}.g.vcf.gz", sample=SAMPLES),
        expand(f"{OUT}/SomVar_VCF/{{sample}}.filtered.vcf.gz", sample=SAMPLES),
        expand(f"{OUT}/bam/flagstat/{{sample}}.txt", sample=SAMPLES)
