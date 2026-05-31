# Code index

This directory is organized by the main MiTRI manuscript workflow. The numbering follows the analysis stages used in the internal `scriptsV4` work directory, with only the cleaned public-facing code retained here.

## Stage map

| Stage | Directory | Purpose |
|---|---|---|
| 00 | `00_preprocessing_and_reads_ratio/` | Reference preparation, metadata preparation, and non-host read ratios. |
| 03 | `03_coverage/` | Coverage generation from microbial alignment outputs. |
| 04 | `04_mitri_activity/` | MiTRI activity and tumour-specific MiTRI calculation. |
| 05 | `05_mtr_region_gene_read_sample/` | MTR vector, region, gene, read, and sample integration. |
| 06 | `06_mtr_peptides/` | MTR peptide and immunogenicity analysis. |
| 07 | `07_peptide_figures/` | Peptide, gene, and species sharing figures. |
| 08 | `08_immunopeptidomics_ms/` | Immunopeptidomics and MS branch analyses. |
| 10 | `10_selective_enrichment/` | Selective enrichment tests for MTR reads. |
| 11 | `11_correlation_and_reads_summary/` | MiTRI read summaries and correlation analyses. |
| 12 | `12_functional_enrichment/` | UniProt-to-GO mapping and GO enrichment. |
| 13 | `13_elispot_selection_and_abundance/` | ELISpot peptide selection and special-cohort abundance figures. |
| 14 | `14_saturation/` | Saturation and accumulation analyses. |
| 15 | `15_single_cell_tcr/` | TCR enrichment statistics. |
| 16 | `16_elispot_response/` | ELISpot response statistics and plots. |
| 99 | `99_manuscript_summaries/` | Cross-stage manuscript summary notebooks. |

## Running modules

Each script or notebook is intended to be run after editing the configuration block near the top of the file. Replace `/path/to/project/` with your local project directory and point the input variables to the matching processed data.

The notebooks are analysis notebooks rather than package APIs. They preserve the manuscript analysis order and are best run section by section.
