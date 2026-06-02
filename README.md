# Microbial transcriptional reprogramming defines tumour-restricted antigens from taxonomically overlapping microbiomes

## Overview

MiTRI (Microbial Transcriptional Reprogramming Index) is a workflow for quantifying microbial transcriptional reprogramming in tumour microenvironments and tracing tumour-restricted microbial transcription to candidate antigens. This repository contains cleaned companion code for the MiTRI manuscript, including MiTRI calculation, MTR region extraction, read and sample integration, peptide analysis, enrichment tests, immunopeptidomics, ELISpot analysis, single-cell TCR summaries, and manuscript figure summaries.

The code is provided for review and reuse. It is not packaged as a one-command pipeline because several steps depend on controlled-access sequencing, local reference indices, and intermediate tables generated on the analysis server.

## Manuscript Status

The MiTRI manuscript is under review and has not been published yet. For that reason, the `paper/` directory is currently a placeholder rather than a final manuscript archive. The manuscript text, figures, and analysis code may still change during peer review. This repository will be updated after publication or when a release-ready manuscript version is available.

## Inputs

The repository starts from processed inputs, including:

- PathSeq or equivalent microbial alignment outputs.
- Per-sample flagstat and microbial read-ratio files.
- Microbial reference genomes and annotation files.
- Coverage bedGraph files and BAM text exports.
- MTR read and sample summary tables.
- HLA binding, pHLA, immunogenicity, MS, ELISpot, and TCR-derived intermediate tables.

Raw sequencing quality control, host filtering, microbial classification, microbial protein discovery, peptide generation, HLA-binding prediction, and host-side preprocessing are maintained in the upstream workflow repository [MimicNeoAI](https://github.com/BioStaCs-public/MimicNeoAI).

All local paths have been anonymized. Replace `/path/to/project/` with the project root on your machine or cluster before running a module.

## Workflow Layout

The `code/` directory follows the main analysis path used in the manuscript.

```text
code/
  00_preprocessing_and_reads_ratio/
      Reference preparation, sample metadata preparation, and non-host read ratio summaries.

  03_coverage/
      Coverage generation helpers for pan-cancer, control, and special cohorts.

  04_mitri_activity/
      Core MiTRI activity calculation and tumour-specific MiTRI summaries.

  05_mtr_region_gene_read_sample/
      MTR vector extraction, MTR region FASTA export, region-to-gene mapping,
      read-level MTR assignment, and sample-level integration.

  06_mtr_peptides/
      MTR peptide filtering, peptide immunogenicity merging, and MiTRI read FASTA export.

  07_peptide_figures/
      Peptide, gene, and species sharing figures and peptide count summaries.

  08_immunopeptidomics_ms/
      Immunopeptidomics and MS analysis notebooks.

  10_selective_enrichment/
      MTR selective enrichment tests and related figures.

  11_correlation_and_reads_summary/
      MiTRI read summaries and MiTRI-immunogenicity correlation analyses.

  12_functional_enrichment/
      UniProt-to-GO mapping and per-cancer GO enrichment.

  13_elispot_selection_and_abundance/
      ELISpot peptide selection and CRC special-cohort abundance figures.

  14_saturation/
      Saturation and accumulation analyses.

  15_single_cell_tcr/
      Single-cell TCR enrichment summaries.

  16_elispot_response/
      ELISpot response statistics and plotting helpers.

  99_manuscript_summaries/
      Figure 1 and Figure 3 summary notebooks that combine outputs from multiple stages.
```

`FiguresOfPaper/` contains lightweight example figure files and small compressed summary inputs retained for manuscript figure reproduction checks.

## Suggested Run Order

Run only the modules needed for the analysis branch you are reproducing.

1. Prepare references and read ratios:
   - `code/00_preprocessing_and_reads_ratio/prepare_reference.ipynb`
   - `code/00_preprocessing_and_reads_ratio/pan_reads_ratio.sh`
   - `code/00_preprocessing_and_reads_ratio/special_reads_ratio.sh`

2. Generate coverage inputs:
   - `code/03_coverage/coverage_pan_cancer.sh`
   - `code/03_coverage/coverage_control.sh`
   - `code/03_coverage/coverage_control_bam2fastq.sh`

3. Calculate MiTRI activity:
   - `code/04_mitri_activity/compute_MiTRI.R`
   - `code/04_mitri_activity/tumour_specific_MiTRI.R`

4. Build MTR regions, genes, reads, and sample summaries:
   - `code/05_mtr_region_gene_read_sample/mtr_region.R`
   - `code/05_mtr_region_gene_read_sample/mtr_genes.ipynb`
   - `code/05_mtr_region_gene_read_sample/mtr_region_reads.R`
   - `code/05_mtr_region_gene_read_sample/mtr_region_per_sample.R`

5. Run peptide, enrichment, MS, GO, ELISpot, TCR, and figure modules as needed:
   - `code/06_mtr_peptides/`
   - `code/07_peptide_figures/`
   - `code/08_immunopeptidomics_ms/`
   - `code/10_selective_enrichment/`
   - `code/11_correlation_and_reads_summary/`
   - `code/12_functional_enrichment/`
   - `code/13_elispot_selection_and_abundance/`
   - `code/14_saturation/`
   - `code/15_single_cell_tcr/`
   - `code/16_elispot_response/`
   - `code/99_manuscript_summaries/`

## Cleaning Notes

For public release, code comments and user-facing messages were converted to English, local server paths were replaced with `/path/to/project/`, and decorative status symbols were replaced with plain text. Variable names, column names, and computational logic were kept unchanged.

## Data Availability

- Sequencing data: PRJNA1301281 (NCBI SRA, controlled access)
- Mass spectrometry proteomics data: PXD074309 (ProteomeXchange, controlled access)
- Antigen analysis: [Zenodo record 19027794](https://zenodo.org/records/19027794)
- Microbial reference: [Zenodo record 18870694](https://zenodo.org/records/18870694)

Large data files and generated analysis outputs are not included in this repository.

## Citation

If you use MiTRI or the associated analysis workflow, please cite:

> Microbial transcriptional reprogramming defines tumour-restricted antigens from taxonomically overlapping microbiomes. (2026)

## License

This software is provided for academic and non-commercial use only. For commercial use, please contact the corresponding authors.

## Contact

For code support or data access questions, contact Yuanwei Zhang (zyuanwei@cmpt.ac.cn).
