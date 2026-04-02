########################################
# CONFIG
########################################
configfile: "config.yaml"

SAMPLES = [line.strip() for line in open(config["sample_file"])]

OUT = config["outdir"]
FQ = config["fq_dir"]

########################################
# FINAL TARGET
########################################
rule all:
    input:
        # Final VCF
        expand(f"{OUT}/SomVar_VCF/{{sample}}.filtered.vcf.gz", sample=SAMPLES),

        # FASTQC RAW
        expand(f"{OUT}/pretrim_fastQC/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        expand(f"{OUT}/pretrim_fastQC/{{sample}}_R2_fastqc.html", sample=SAMPLES),

        # FASTQC TRIM ✅ FIXED
        expand(f"{OUT}/post_trim_fastQC/{{sample}}_R1_val_1_fastqc.html", sample=SAMPLES),
        expand(f"{OUT}/post_trim_fastQC/{{sample}}_R2_val_2_fastqc.html", sample=SAMPLES),

        # MULTIQC
        f"{OUT}/pretrim_fastQC/multiqc_report.html",
        f"{OUT}/post_trim_fastQC/multiqc_report.html"

########################################
# 1. FASTQC RAW
########################################
rule fastqc_raw:
    input:
        r1 = f"{FQ}/{{sample}}_R1.fastq.gz",
        r2 = f"{FQ}/{{sample}}_R2.fastq.gz"
    output:
        r1 = f"{OUT}/pretrim_fastQC/{{sample}}_R1_fastqc.html",
        r2 = f"{OUT}/pretrim_fastQC/{{sample}}_R2_fastqc.html"
    threads: config["threads"]["fastqc"]
    shell:
        """
        mkdir -p {OUT}/pretrim_fastQC
        fastqc -t {threads} -o {OUT}/pretrim_fastQC {input.r1} {input.r2}
        """

########################################
# MULTIQC RAW
########################################
rule multiqc_raw:
    input:
        expand(f"{OUT}/pretrim_fastQC/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        expand(f"{OUT}/pretrim_fastQC/{{sample}}_R2_fastqc.html", sample=SAMPLES)
    output:
        f"{OUT}/pretrim_fastQC/multiqc_report.html"
    shell:
        """
        multiqc {OUT}/pretrim_fastQC \
        -o {OUT}/pretrim_fastQC \
        -n multiqc_report.html
        """

########################################
# 2. TRIMMING
########################################
rule trim_galore:
    input:
        r1 = f"{FQ}/{{sample}}_R1.fastq.gz",
        r2 = f"{FQ}/{{sample}}_R2.fastq.gz",
        qc = f"{OUT}/pretrim_fastQC/{{sample}}_R1_fastqc.html"
    output:
        r1 = f"{OUT}/TrimGalore/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{OUT}/TrimGalore/{{sample}}_R2_val_2.fq.gz"
    threads: config["threads"]["trim"]
    shell:
        """
        mkdir -p {OUT}/TrimGalore
        trim_galore --paired --cores {threads} \
        -o {OUT}/TrimGalore {input.r1} {input.r2}
        """

########################################
# 3. FASTQC AFTER TRIM
########################################
rule fastqc_trim:
    input:
        r1 = f"{OUT}/TrimGalore/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{OUT}/TrimGalore/{{sample}}_R2_val_2.fq.gz"
    output:
        r1 = f"{OUT}/post_trim_fastQC/{{sample}}_R1_val_1_fastqc.html",
        r2 = f"{OUT}/post_trim_fastQC/{{sample}}_R2_val_2_fastqc.html",
        z1 = f"{OUT}/post_trim_fastQC/{{sample}}_R1_val_1_fastqc.zip",
        z2 = f"{OUT}/post_trim_fastQC/{{sample}}_R2_val_2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    shell:
        """
        mkdir -p {OUT}/post_trim_fastQC
        fastqc -t {threads} -o {OUT}/post_trim_fastQC {input.r1} {input.r2}
        """

########################################
# MULTIQC POST TRIM
########################################
rule multiqc_trim:
    input:
        expand(f"{OUT}/post_trim_fastQC/{{sample}}_R1_val_1_fastqc.html", sample=SAMPLES),
        expand(f"{OUT}/post_trim_fastQC/{{sample}}_R2_val_2_fastqc.html", sample=SAMPLES)
    output:
        f"{OUT}/post_trim_fastQC/multiqc_report.html"
    shell:
        """
        multiqc {OUT}/post_trim_fastQC -o {OUT}/post_trim_fastQC -n multiqc_report.html
        """

########################################
# 4. ALIGNMENT
########################################
rule align:
    input:
        r1 = f"{OUT}/TrimGalore/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{OUT}/TrimGalore/{{sample}}_R2_val_2.fq.gz"
    output:
        f"{OUT}/bam/{{sample}}_sorted.bam"
    threads: config["threads"]["dragmap"]
    params:
        ref = config["ref_genome"],
        samthreads = config["threads"]["samthreads"]
    shell:
        """
        set -o pipefail
        dragen-os -r {params.ref} \
        -1 {input.r1} -2 {input.r2} \
        --num-threads {threads} \
        --RGSM {wildcards.sample} | \
        samtools sort -@ {params.samthreads} -o {output}
        """

########################################
# 5. MARK DUPLICATES
########################################
rule markdup:
    input:
        f"{OUT}/bam/{{sample}}_sorted.bam"
    output:
        bam = f"{OUT}/bam/{{sample}}_sorted_MD.bam",
        metrics = f"{OUT}/picard_info/{{sample}}.txt"
    shell:
        """
        mkdir -p {OUT}/picard_info
        picard MarkDuplicates \
        I={input} O={output.bam} \
        M={output.metrics} \
        ASSUME_SORT_ORDER=coordinate
        """

########################################
# 6. INDEX BAM
########################################
rule index:
    input:
        bam = f"{OUT}/bam/{{sample}}_sorted_MD.bam"
    output:
        bai = f"{OUT}/bam/{{sample}}_sorted_MD.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

########################################
# 7. FLAGSTAT QC
########################################
rule flagstat:
    input:
        f"{OUT}/bam/{{sample}}_sorted_MD.bam"
    output:
        f"{OUT}/bam/flagstat/{{sample}}.txt"
    shell:
        """
        mkdir -p {OUT}/bam/flagstat
        samtools flagstat {input} > {output}
        """

########################################
# 8. MUTECT2
########################################
rule mutect2:
    input:
        bam = f"{OUT}/bam/{{sample}}_sorted_MD.bam",
        bai = f"{OUT}/bam/{{sample}}_sorted_MD.bam.bai"
    output:
        vcf = f"{OUT}/SomVar_VCF/{{sample}}.vcf.gz",
        f1r2 = f"{OUT}/SomVar_VCF/{{sample}}.f1r2.tar.gz"
    threads: config["threads"]["gatk"]
    params:
        ref = config["mut_ref"],
        germ = config["germline"],
        pon = config["pon"]
    shell:
        """
        mkdir -p {OUT}/SomVar_VCF
        gatk Mutect2 \
        -R {params.ref} \
        -I {input.bam} \
        --germline-resource {params.germ} \
        --panel-of-normals {params.pon} \
        --f1r2-tar-gz {output.f1r2} \
        -O {output.vcf}
        """

########################################
# 9. LEARN ORIENTATION
########################################
rule learn_orientation:
    input:
        f"{OUT}/SomVar_VCF/{{sample}}.f1r2.tar.gz"
    output:
        f"{OUT}/SomVar_VCF/{{sample}}.artifact.tar.gz"
    shell:
        """
        gatk LearnReadOrientationModel -I {input} -O {output}
        """

########################################
# 10. FILTER
########################################
rule filter:
    input:
        vcf = f"{OUT}/SomVar_VCF/{{sample}}.vcf.gz",
        art = f"{OUT}/SomVar_VCF/{{sample}}.artifact.tar.gz"
    output:
        f"{OUT}/SomVar_VCF/{{sample}}.filtered.vcf.gz"
    params:
        ref = config["mut_ref"]
    shell:
        """
        gatk FilterMutectCalls \
        -V {input.vcf} \
        --orientation-bias-artifact-priors {input.art} \
        -R {params.ref} \
        -O {output}
        """
