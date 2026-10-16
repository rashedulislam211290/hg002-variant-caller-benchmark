# Benchmarking BCFtools and FreeBayes against the GIAB HG002 truth set

This repository contains the workflow used to compare two open-source germline
small-variant callers, **BCFtools** and **FreeBayes**, on identical input data,
and to score both against the **Genome in a Bottle (GIAB) HG002 v4.2.1**
benchmark with **hap.py**, following the GA4GH best-practice framework.

The analysis is restricted to **chromosome 20 (GRCh38)** so that the full workflow
runs in under an hour on a single compute node.

## Workflow

| Step | Script | Where it runs | What it does |
|---|---|---|---|
| 1 | `scripts/01_download_reference.sh` | login node (needs internet) | Downloads the GRCh38 no-alt analysis set and indexes it with `samtools faidx` |
| 2 | `scripts/02_download_data.sh` | login node (needs internet) | Downloads the HG002 chr20 BAM and the GIAB v4.2.1 truth VCF and BED |
| 3–5 | `scripts/03_call_and_benchmark.pbs` | compute node (PBS job) | Calls variants with BCFtools and FreeBayes, normalises with `bcftools norm`, benchmarks with hap.py |
| 6 | `scripts/04_make_figure.py` | anywhere with Python | Builds Figure 1 and a performance table from the hap.py output |

## Data

All inputs are public and are downloaded by the scripts.

| Input | Source |
|---|---|
| Reference: GRCh38 no-alt analysis set (`GCA_000001405.15`) | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/) |
| Reads: HG002, Illumina NovaSeq PCR-free, ~35×, aligned to GRCh38 no-alt, chr20 | [DeepVariant case-study data](https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam) |
| Truth set: GIAB HG002 v4.2.1 benchmark VCF + benchmark regions BED | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/) |

## Requirements

- Linux with conda (Miniconda or Mambaforge)
- A PBS Professional cluster for step 3–5 (any Linux machine works if you run the script with `bash` instead of `qsub`)
- About 15 GB of free disk space and 16 GB of RAM

Create the environments once:

```bash
conda env create -f envs/ngs.yml     # samtools, bcftools, htslib, freebayes
conda env create -f envs/happy.yml   # hap.py (Python 2.7, kept separate)
conda env create -f envs/plot.yml    # pandas, matplotlib for the figure
```

## Usage

Run everything from the repository root:

```bash
conda activate ngs
bash scripts/01_download_reference.sh
bash scripts/02_download_data.sh

qsub scripts/03_call_and_benchmark.pbs      # check with: qstat -u $USER

conda activate plot
python scripts/04_make_figure.py --results results --out figures
```

The queue name (`thin`) and resource requests in the PBS header may need to be
changed for your cluster.

## Methods summary

- **BCFtools:** `bcftools mpileup -a AD,DP | bcftools call -mv` (multiallelic model, variant sites only)
- **FreeBayes:** default parameters
- No post-calling filters were applied to either call set.
- **Normalisation:** `bcftools norm -f <ref> -m -any` (split multiallelic records, left-align indels)
- **Benchmarking:** `hap.py` (xcmp engine) restricted to the GIAB benchmark regions and to chr20 (`-l chr20`)

## Results

Performance within the GIAB v4.2.1 benchmark regions on chr20 (71,333 SNVs and 11,256 indels):

| Variant type | Caller | TP | FN | FP | Recall | Precision | F1 |
|---|---|---|---|---|---|---|---|
| SNV | BCFtools | 70,894 | 439 | 790 | 0.9938 | **0.9890** | **0.9914** |
| SNV | FreeBayes | 70,991 | 342 | 2,341 | **0.9952** | 0.9681 | 0.9815 |
| Indel | BCFtools | 10,552 | 704 | 596 | 0.9375 | 0.9464 | 0.9419 |
| Indel | FreeBayes | 10,887 | 369 | 352 | **0.9672** | **0.9687** | **0.9679** |

BCFtools was more accurate for SNVs, mainly because it produced fewer false
positives; FreeBayes was more accurate for indels, with both higher recall and
higher precision.

![Figure 1](figures/Figure1.png)

## Repository layout

```
envs/        conda environment files
scripts/     workflow scripts (steps 1-6)
results/     hap.py summary files (*.summary.csv)
figures/     Figure 1
```

Large data files (FASTA, BAM, VCF) are not tracked; see `.gitignore`.

## References

1. Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience* 2021;10:giab008.
2. Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv* 2012;1207.3907.
3. Wagner J, et al. Benchmarking challenging small variants with linked and long reads. *Cell Genomics* 2022;2:100128.
4. Krusche P, et al. Best practices for benchmarking germline small-variant calls in human genomes. *Nat Biotechnol* 2019;37:555–560.

## License

MIT, see [LICENSE](LICENSE).
