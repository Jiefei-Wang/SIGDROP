# SIGDROP: candidate PSA declines

`sigdrop.R` provides the standalone base-R function `find_sigdrop()` for the paper "Developing an OMOP-Standardized Prostate Cancer Database and Improving Data Quality Using NLP and PSA-Based Algorithms".

## Original paper

Adapted from Algorithm 1 (SIGDROP) in:

Batouche AO, Czeizler E, Koskinen M, Mirtti T, Rannikko AS. (2023). *Synergizing Data Imputation and Electronic Health Records for Advancing Prostate Cancer Research: Challenges, and Practical Applications.* arXiv:2311.02086. [Paper and original pseudocode](https://arxiv.org/html/2311.02086v1).

This is the study's R adaptation, not the original authors' source code.

## Usage

From the repository root:

```r
source("sigdrop.R")
find_sigdrop(psa = c(8, 0.05), days = c(0, 30))
#      drop_day nadir_day PSAmax PSAmin
# item        0        30      8   0.05
```

## Modifications and implementation details

- **Multiple candidates:** Save each first threshold-crossing decline and reset the peak to the following measurement (`j + 2`). The published pseudocode returns one result after scanning the series.
- **Division-by-zero protection:** Calculate the relative decline using `max(PSAmax, 0.0001)` as the denominator.
- **Time and boundary handling:** Express the 12-month reset interval as 365 days, calculate elapsed time from the current next measurement, and check bounds before accessing `j + 2`. Equal peaks also reset the peak date (`>=`, versus `>` in the original pseudocode).
- **Input/output handling:** Sort paired measurements by time, require equal vector lengths, return `NULL` for fewer than two measurements, and include the peak value in each output row.

The significant-decline thresholds are unchanged: a relative decline of at least 75% with an absolute decline of at least 3 ng/mL, or at least 50% with an absolute decline of at least 4 ng/mL.
