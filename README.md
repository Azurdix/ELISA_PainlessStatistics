# ELISA Pipeline

A compact R workflow for analysing **ELISA-based protein concentration data**.

This script was developed as part of my PhD project and is designed to work with an external R file containing the input data vectors.

## Overview

The script:

- loads ELISA data from a separate R script,
- converts all values to a common concentration scale (`µg/mL`),
- reshapes the data into long format,
- compares `EGGS` and `SILK` samples against the `BLANK`,
- applies multiple-testing correction,
- generates summary tables,
- creates publication-style plots,
- adds significance labels,
- saves final outputs as high-resolution files.

## Input

The script requires an external R file, for example:

`source(".../protein_matrices.R")`

That file should contain:

- BLANK
- protein-specific vectors such as:
- PRDE_LIZO
- PRDE_DEFE
- PRDE_CECRO
- PRDS_LIZO
- PRDS_DEFE
- PRDS_CECRO
- PTSDE_LIZO
- PTSDE_DEFE
- PTSDE_CECRO
- PTSDS_LIZO
- PTSDS_DEFE
- PTSDS_CECRO

## Default groups

Sample types:
- EGGS
- SILK

Spider species:
- PRD - _Pardosa lugubris_
- PTSD - _Parasteatoda tepidariorum_

Protein groups:
- LIZO – Lysozyme-like peptides
- DEFE – Defensin-like peptides
- CECRO – Cecropin-like peptides

## Statistics
By default, the script uses:

- one-sided Wilcoxon rank-sum test
- alternative hypothesis: group > blank
- Holm correction for multiple comparisons

Comparisons are performed separately within each protein × spider species combination.

## Output

The script produces:

- a formatted statistical summary table,
- a CSV file with test results,
- a raw-value plot with blank reference line and IQR band,
- a blank-corrected plot,
- significance annotations for both figures.

Saved files:

- protein_statistics_vs_blank_one_sided.csv
- protein_concentration_raw_with_blank_line.tiff
- protein_concentration_blank_corrected.tiff

## Notes
This script was written to be flexible and easy to adapt to expanding ELISA datasets.
Some internal comments may still contain Polish.
