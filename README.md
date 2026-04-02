# variant_analysis_pipeline
A modular Snakemake workflow implementing GATK best practices for germline and somatic variant calling, featuring HaplotypeCaller, Mutect2, FastQC, Trim Galore, DRAGEN alignment, and Picard MarkDuplicates.

# 🧬 Reproducible Variant Calling Pipelines (Snakemake)

This repository provides **modular Snakemake pipelines** for:

- 🔬 Somatic Variant Calling (Mutect2)
- 🧬 Germline Variant Calling (HaplotypeCaller + Joint Genotyping)

---

## 📂 File Structure

### 🔬 Somatic Variant Pipeline
- `somatic_variant.smk`
- `somatic_config.yaml`
- `requirement.txt`

### 🧬 Germline Variant Pipeline
- `germline_variant.smk`
- `germline_config.yaml`
- `requirement.txt`

---

## ⚙️ Installation

```bash
conda create -n variant_env --file requirement.txt
conda activate variant_env

## ▶️ Running the Pipeline

```bash
snakemake --cores 6 --configfile <config.yaml> --rerun-incomplete --keep-going --latency-wait 20
