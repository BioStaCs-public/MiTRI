# Quantifying microbial transcriptional reprogramming enables tumour-restricted peptide discovery

## Overview

MiTRI (Microbial Transcriptional Reprogramming Index) is a computational framework for quantifying microbial transcriptional reprogramming in tumour microenvironments. This repository contains the cleaned analysis code used for MiTRI quantification, downstream enrichment, immunopeptidomics, ELISpot analysis, single-cell/TCR summaries, saturation analysis, and manuscript figure generation.

## Scope and Inputs

This repository starts from processed microbial alignment, read-count, peptide-prediction, and annotation outputs. Raw sequencing quality control, host/vector filtering, microbial taxonomic classification, microbial protein or peptide discovery, HLA-binding prediction, and host genetics/transcriptomics preprocessing are maintained in the upstream pipeline repository [MimicNeoAI](https://github.com/BioStaCs-public/MimicNeoAI).

Please run MimicNeoAI, or prepare equivalent processed inputs, before using the MiTRI modules here.

## Repository Structure

```text
code/
  01_reads_ratio/              Reference preparation and read-ratio summaries
  02_coverage/                 Coverage calculation helpers
  03_mitri_calculation/        Core MiTRI score calculation and read statistics
  04_mtr_region/               MTR region and MTR gene quantification
  05_selective_enrichment/     Selective transcriptional enrichment analysis
  06_mtr_peptides/             MTR peptide extraction, filtering, and FASTA export
  07_functional_enrichment/    GO enrichment analysis
  08_immunopeptidomics/        Immunopeptidomics and UniProt-GO mapping analyses
  09_single_cell_tcr/          Single-cell TCR enrichment summaries
  10_elispot/                  ELISpot response analysis
  11_correlation/              MiTRI and immunogenicity correlation analysis
  12_saturation/               Saturation and accumulation analyses
  13_figures/                  Figure-generation scripts and summary notebooks
FiguresOfPaper/                Lightweight figure examples and plotted outputs
LICENSE
README.md
```

## Data Availability

- Sequencing data: PRJNA1301281 (NCBI SRA, controlled access)
- Mass spectrometry proteomics data: PXD074309 (ProteomeXchange, controlled access)
- Antigen analysis: [Zenodo record 19027794](https://zenodo.org/records/19027794)
- Microbial reference: [Zenodo record 18870694](https://zenodo.org/records/18870694)

## Citation

If you use MiTRI or the associated analysis workflow, please cite:

> Quantifying microbial transcriptional reprogramming enables tumour-restricted peptide discovery. (2026)

## License

This software is provided for academic and non-commercial use only. For commercial use, please contact the corresponding authors.

## Contact

For code support or data access questions, contact Yuanwei Zhang (zyuanwei@cmpt.ac.cn).

## Acknowledgments

We gratefully acknowledge the Hefei Advanced Computing Center for computational resources. Funding support is described in the manuscript.