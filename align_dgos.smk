# ============================================================
# align_dgos.smk — DRAGEN-OS alignment + QC pipeline
# Converted from the reference bash script.
#
# Run with, e.g.:
#   snakemake -s align_dgos.smk --configfile config.yaml -c 32 -p
#
# Notes on the conversion:
#   * Conda/env activation is NOT handled here — activate your
#     environment yourself before running `snakemake`.
#   * nohup/background execution and manual `parallel -j N` fan-out are
#     no longer needed — Snakemake parallelizes automatically across
#     samples/rules based on the `-c/--cores` you give it and each
#     rule's `threads:` value.
#   * `dragen-os | samtools sort` is kept as a single piped shell
#     command inside one rule (as in the original), so the combined
#     rule's thread request is align + sort threads.
# ============================================================

import os

configfile: "config.yaml"

# ----------------- Samples -----------------
with open(config["samples_file"]) as fh:
    SAMPLES = [line.strip() for line in fh if line.strip()]

# ----------------- Paths -----------------
FQ      = config["fq_data"]
RESULT  = config["result"]
REF     = config["ref_genome"]

PRETRIM_QC  = f"{RESULT}/pretrim_fastQC"
TRIM        = f"{RESULT}/TrimGalore"
POSTTRIM_QC = f"{RESULT}/post_trim_fastQC"
BAM         = f"{RESULT}/bam"
PICARD      = f"{RESULT}/picard_info"
VCF         = f"{RESULT}/SomVar_VCF"

BAM_FASTQC = f"{BAM}/fastqc"
QUALIMAP   = f"{BAM}/qualimap"
FLAGSTAT   = f"{BAM}/flagstat"
STATS      = f"{BAM}/stats"

T = config["threads"]
TG = config["trim_galore"]

wildcard_constraints:
    sample = "|".join([s.replace(".", r"\.") for s in SAMPLES]) if SAMPLES else ".*"


# ============================================================
# rule all — final targets
# ============================================================
rule all:
    input:
        f"{PRETRIM_QC}/multiqc_report.html",
        f"{POSTTRIM_QC}/multiqc_report.html",
        expand(f"{BAM}/{{sample}}_dragen_sorted_MD.bam.bai", sample=SAMPLES),
        expand(f"{FLAGSTAT}/{{sample}}.flagstat.txt", sample=SAMPLES),
        expand(f"{STATS}/{{sample}}.stats.txt", sample=SAMPLES),
        expand(f"{QUALIMAP}/{{sample}}_qualimap", sample=SAMPLES),
        f"{RESULT}/bam_qc.html",
        f"{BAM}/multiqc_report.html",


# ============================================================
# 1. FastQC: pre-trim quality check
# ============================================================
rule fastqc_pretrim:
    input:
        r1 = f"{FQ}/{{sample}}_R1.fastq.gz",
        r2 = f"{FQ}/{{sample}}_R2.fastq.gz",
    output:
        html1 = f"{PRETRIM_QC}/{{sample}}_R1_fastqc.html",
        zip1  = f"{PRETRIM_QC}/{{sample}}_R1_fastqc.zip",
        html2 = f"{PRETRIM_QC}/{{sample}}_R2_fastqc.html",
        zip2  = f"{PRETRIM_QC}/{{sample}}_R2_fastqc.zip",
    threads: T["fastqc"]
    shell:
        "fastqc --outdir {PRETRIM_QC} -t {threads} {input.r1} {input.r2}"


rule multiqc_pretrim:
    input:
        expand(f"{PRETRIM_QC}/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        expand(f"{PRETRIM_QC}/{{sample}}_R2_fastqc.zip", sample=SAMPLES),
    output:
        f"{PRETRIM_QC}/multiqc_report.html",
    shell:
        "multiqc -o {PRETRIM_QC} {PRETRIM_QC}"


# ============================================================
# 2. Trim Galore
# ============================================================
rule trim_galore:
    input:
        r1 = f"{FQ}/{{sample}}_R1.fastq.gz",
        r2 = f"{FQ}/{{sample}}_R2.fastq.gz",
    output:
        r1 = f"{TRIM}/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{TRIM}/{{sample}}_R2_val_2.fq.gz",
        report1 = f"{TRIM}/{{sample}}_R1.fastq.gz_trimming_report.txt",
        report2 = f"{TRIM}/{{sample}}_R2.fastq.gz_trimming_report.txt",
    threads: T["trim_galore"]
    shell:
        "trim_galore --paired --cores {threads} --quality {TG[quality]} "
        "--fastqc --phred33 --length {TG[length]} --stringency {TG[stringency]} "
        "-o {TRIM} {input.r1} {input.r2}"


# ============================================================
# 3. FastQC: post-trim quality check
# ============================================================
rule fastqc_posttrim:
    input:
        r1 = f"{TRIM}/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{TRIM}/{{sample}}_R2_val_2.fq.gz",
    output:
        html1 = f"{POSTTRIM_QC}/{{sample}}_R1_val_1_fastqc.html",
        zip1  = f"{POSTTRIM_QC}/{{sample}}_R1_val_1_fastqc.zip",
        html2 = f"{POSTTRIM_QC}/{{sample}}_R2_val_2_fastqc.html",
        zip2  = f"{POSTTRIM_QC}/{{sample}}_R2_val_2_fastqc.zip",
    threads: T["fastqc"]
    shell:
        "fastqc --outdir {POSTTRIM_QC} -t {threads} {input.r1} {input.r2}"


rule multiqc_posttrim:
    input:
        expand(f"{POSTTRIM_QC}/{{sample}}_R1_val_1_fastqc.zip", sample=SAMPLES),
        expand(f"{POSTTRIM_QC}/{{sample}}_R2_val_2_fastqc.zip", sample=SAMPLES),
    output:
        f"{POSTTRIM_QC}/multiqc_report.html",
    shell:
        "multiqc -o {POSTTRIM_QC} {POSTTRIM_QC}"


# ============================================================
# 4. Alignment: dragen-os | samtools sort
# ============================================================
rule align:
    input:
        r1 = f"{TRIM}/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{TRIM}/{{sample}}_R2_val_2.fq.gz",
    output:
        bam = f"{BAM}/{{sample}}_dragen_sorted.bam",
    threads: T["align"] + T["sort"]
    params:
        ref = REF,
        align_threads = T["align"],
        sort_threads = T["sort"],
        sort_mem = config["mem"]["sort_per_thread"],
    shell:
        "dragen-os -r {params.ref} "
        "-1 {input.r1} -2 {input.r2} "
        "--num-threads {params.align_threads} "
        "--RGID illumina "
        "--Aligner.smith-waterman-method dragen "
        "--RGSM {wildcards.sample} | "
        "samtools sort -@ {params.sort_threads} -m {params.sort_mem} -o {output.bam}"


# ============================================================
# 5. MarkDuplicates
# ============================================================
rule mark_duplicates:
    input:
        bam = f"{BAM}/{{sample}}_dragen_sorted.bam",
    output:
        bam = f"{BAM}/{{sample}}_dragen_sorted_MD.bam",
        metrics = f"{PICARD}/{{sample}}_MarkDup_metrics.txt",
    threads: T["markdup"]
    shell:
        "picard MarkDuplicates "
        "-I {input.bam} "
        "-O {output.bam} "
        "-M {output.metrics} "
        "--REMOVE_DUPLICATES false "
        "--ASSUME_SORT_ORDER coordinate"


# ============================================================
# 6. Index
# ============================================================
rule index_bam:
    input:
        bam = f"{BAM}/{{sample}}_dragen_sorted_MD.bam",
    output:
        bai = f"{BAM}/{{sample}}_dragen_sorted_MD.bam.bai",
    threads: T["index"]
    shell:
        "samtools index -@ {threads} {input.bam}"


# ============================================================
# 7. QC post-alignment: FastQC, Qualimap, flagstat, stats
# ============================================================
rule fastqc_bam:
    input:
        bam = f"{BAM}/{{sample}}_dragen_sorted_MD.bam",
        bai = f"{BAM}/{{sample}}_dragen_sorted_MD.bam.bai",
    output:
        html = f"{BAM_FASTQC}/{{sample}}_dragen_sorted_MD_fastqc.html",
        zip  = f"{BAM_FASTQC}/{{sample}}_dragen_sorted_MD_fastqc.zip",
    threads: T["qc"]
    shell:
        "fastqc -t {threads} -o {BAM_FASTQC} {input.bam}"


rule qualimap_bam:
    input:
        bam = f"{BAM}/{{sample}}_dragen_sorted_MD.bam",
        bai = f"{BAM}/{{sample}}_dragen_sorted_MD.bam.bai",
    output:
        outdir = directory(f"{QUALIMAP}/{{sample}}_qualimap"),
    threads: T["qc"]
    shell:
        "qualimap bamqc -bam {input.bam} -outdir {output.outdir}"


rule flagstat_bam:
    input:
        bam = f"{BAM}/{{sample}}_dragen_sorted_MD.bam",
        bai = f"{BAM}/{{sample}}_dragen_sorted_MD.bam.bai",
    output:
        f"{FLAGSTAT}/{{sample}}.flagstat.txt",
    shell:
        "samtools flagstat {input.bam} > {output}"


rule stats_bam:
    input:
        bam = f"{BAM}/{{sample}}_dragen_sorted_MD.bam",
        bai = f"{BAM}/{{sample}}_dragen_sorted_MD.bam.bai",
    output:
        f"{STATS}/{{sample}}.stats.txt",
    shell:
        "samtools stats {input.bam} > {output}"


# ============================================================
# 8. MultiQC over the whole bam/ directory (two reports, as in original)
# ============================================================
rule multiqc_bam_qc:
    input:
        expand(f"{BAM_FASTQC}/{{sample}}_dragen_sorted_MD_fastqc.zip", sample=SAMPLES),
        expand(f"{QUALIMAP}/{{sample}}_qualimap", sample=SAMPLES),
        expand(f"{FLAGSTAT}/{{sample}}.flagstat.txt", sample=SAMPLES),
        expand(f"{STATS}/{{sample}}.stats.txt", sample=SAMPLES),
    output:
        f"{RESULT}/bam_qc.html",
    shell:
        "multiqc {BAM} -o {RESULT} -n bam_qc.html"


rule multiqc_bam:
    input:
        f"{RESULT}/bam_qc.html",
    output:
        f"{BAM}/multiqc_report.html",
    shell:
        "multiqc {BAM} -o {BAM}"
