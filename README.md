# Variant Analysis Pipeline
A modular Snakemake workflow implementing GATK best practices for germline and somatic variant calling, featuring HaplotypeCaller, Mutect2, FastQC, Trim Galore, DRAGEN alignment, and Picard MarkDuplicates.

#  Reproducible Variant Calling Pipelines (Snakemake)

This repository provides **modular Snakemake pipelines** for:

-  Somatic Variant Calling (Mutect2)
-  Germline Variant Calling (HaplotypeCaller + Joint Genotyping)

---

##  File Structure

###  Somatic Variant Pipeline
- `somatic_variant.smk`
- `somatic_config.yaml`
- `requirement.txt`

###  Germline Variant Pipeline
- `germline_variant.smk`
- `germline_config.yaml`
- `requirement.txt`

---
---

## ⚙️ Configuration File

The pipeline is controlled using a YAML configuration file  
(e.g. `somatic_variant_config.yaml` or `germline_variant_config.yaml`).

### 📌 Example Configuration

```yaml
fq_dir: "/path/to/raw_fastq"
outdir: "/path/to/output_directory"

sample_file: "/path/to/sample_list.txt"

ref_genome: "/path/to/dragmap_reference/"
mut_ref: "/path/to/hg38.fa"
germline: "/path/to/germline_resource.vcf.gz"
pon: "/path/to/panel_of_normals.vcf.gz"

threads:
  fastqc: 4
  trim: 4
  dragmap: 4
  samthreads: 2
  gatk: 4
```

---

### 🧾 Parameter Description

| Parameter | Description |
|----------|------------|
| `fq_dir` | Directory containing input FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`) |
| `outdir` | Output directory where all results will be stored |
| `sample_file` | Text file with sample names (one per line) |
| `ref_genome` | DRAGEN reference genome directory (used for alignment) |
| `mut_ref` | Reference FASTA file (used by GATK tools) |
| `germline` | Germline resource VCF (e.g. gnomAD) for Mutect2 |
| `pon` | Panel of Normals (PoN) VCF file |
| `threads.fastqc` | Threads for FastQC |
| `threads.trim` | Threads for Trim Galore |
| `threads.dragmap` | Threads for DRAGEN alignment |
| `threads.samthreads` | Threads for samtools sorting |
| `threads.gatk` | Threads for GATK tools |

---

### 📄 Sample File Format

The `sample_file` should contain one sample per line:

```txt
Sample1
Sample2
Sample3
```

> ⚠️ Make sure FASTQ files follow this naming pattern:
> `Sample1_R1.fastq.gz` and `Sample1_R2.fastq.gz`

---

### ⚠️ Notes

- Paths must be **absolute paths** (recommended for cluster/HPC runs)
- Ensure all reference files are **indexed**:
  - `.fai` for FASTA
  - `.dict` for GATK
  - `.bai` for BAM
- Thread settings can be adjusted based on your system resources

---

## ⚙️ Installation

```bash
conda create -n variant_env --file requirement.txt
conda activate variant_env
```

---

##  Running the Pipeline

```bash
snakemake -s <somatic_variant.smk OR germline_variant.smk> \
  --cores 6 \
  --configfile <somatic_variant_config.yaml OR germline_variant_config.yaml> \
  --rerun-incomplete \
  --keep-going \
  --latency-wait <#>
```

---

##   Example

```bash
snakemake -s somatic_variant.smk \
  --cores 6 \
  --configfile somatic_variant_config.yaml \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 20
```

### 🌐 Reference Resources

Pre-built **Panel of Normals (PoN)** and **germline resource files** are available from GATK Best Practices:

👉 https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38

---
