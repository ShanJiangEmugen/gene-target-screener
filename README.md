# U7 Target Screener

A small Biopython-based tool for scanning U7 antisense target sequences
across species or between different transcript regions.

This repo contains the code I use in my AAV–U7 splicing correction
projects to identify conserved, low-repeat subsequences that are
compatible with U7-snRNA–mediated splice modulation.

## Features

- Sliding-window extraction of candidate subsequences (default 18–35 bp)
- Filters out CTGCTGCTG-rich regions (optional to extend later)
- Pairwise alignment using Biopython `PairwiseAligner` with BLASTN matrix
- Simple conservation and score-based filtering
- Batch mode: screen many pairwise comparisons from a single CSV file
- Outputs results as CSV files for downstream cloning/construct design

## Installation

```bash
git clone https://github.com/<your-username>/u7-target-screener.git
cd u7-target-screener
pip install -r requirements.txt

