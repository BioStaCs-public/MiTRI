# Quantifying microbial transcriptional reprogramming enables tumour-restricted peptide discovery

## Overview
**MiTRI** (Microbial Transcriptional Reprogramming Index) is a specialized computational framework developed to quantify microbial transcriptional reprogramming in tumor microenvironments. This repository contains the core analytical scripts used to calculate MiTRI values, perform selective transcriptional enrichment analysis, and reproduce the figures presented in our study.

⚠️ **Pipeline Architecture & Prerequisites:**
This repository (`MiTRI`) focuses on the specific quantification of transcriptional reprogramming and downstream visualizations. 
The comprehensive upstream processing—including raw data quality control, host/vector filtering, microbial taxonomic classification, microbial protein/peptide discovery, HLA-binding prediction, and host genetics/transcriptomics—is hosted in our foundational pipeline repository: **[MimicNeoAI](https://github.com/BioStaCs-public/MimicNeoAI)**. 

Please ensure you have processed your raw sequencing datasets through the `MimicNeoAI` pipeline (or have equivalent processed data) before running the `MiTRI` modules.

### Key Features Included in this Repository
- **Prevalence Filtering:** Filter microbial taxa to retain only those with reliable detection rates (e.g., ≥ 50% prevalence) across tumor/peritumor sample cohorts.
- **Microbial Transcriptional Reprogramming (MiTRI) Quantification:** - Precise target bacteria re-alignment and stringent identity filtering.
  - Genome-wide transcriptional profiling and background noise normalization.
  - Calculation of the MiTRI score to quantify the magnitude of tumor-restricted transcriptional divergence.
- **Selective Transcriptional Enrichment Analysis:** Statistical evaluation to confirm that observed reprogramming is driven by selective transcriptional regulation rather than passive genomic abundance.
- **Figure Plotting & Visualization:** R/Python scripts and necessary processed data required to reproduce the main figures and extended data figures from the manuscript.

## Citation
If you use MiTRI or our methodology in your research, please cite our paper:
> "Quantifying microbial transcriptional reprogramming enables tumour-restricted peptide discovery." (2026)

## Data Availability
- **Sequencing Data:** PRJNA1301281 (NCBI SRA, controlled access)
- **Mass Spectrometry Proteomics Data:** PXD074309 (ProteomeXchange, controlled access)
- **Antigen Analysis:** [Zenodo 10.5281/zenodo.15041208](https://zenodo.org/record/15041208)
- **Microbial Reference:** [Zenodo 10.5281/zenodo.15043056](https://zenodo.org/record/15043056) 
  
## License
This software is provided for academic and non-commercial use only. 
For commercial use, please contact the corresponding authors.

## Contact
For questions, code-related support, or data access requests, please contact the corresponding authors: 
- Yuanwei Zhang (zyuanwei@cmpt.ac.cn)

## Acknowledgments
We gratefully acknowledge the Hefei Advanced Computing Center for providing computational resources. This work was supported by the National Natural Science Foundation for Young Scholars of China (82171601), the Anhui Provincial Key Research and Development Plan (2022e07020054), and other funding agencies as detailed in the manuscript.
