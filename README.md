<h1 align="left">Gene Target Screener</h1>

<p align="left">
  <img src="https://img.shields.io/badge/status-active-brightgreen" alt="Project Status">
  <img src="https://img.shields.io/badge/version-0.1.0-blue" alt="Version">
  <img src="https://img.shields.io/github/license/ShanJiangEmugen/gene-target-screener" alt="License">
  <img src="https://img.shields.io/github/stars/ShanJiangEmugen/gene-target-screener?style=social" alt="GitHub Stars">
</p>

<p align="left">
  <img src="https://img.shields.io/badge/python-3.8%2B-blue" alt="Python Version">
  <img src="https://img.shields.io/github/last-commit/ShanJiangEmugen/gene-target-screener" alt="Last Commit">
  <img src="https://img.shields.io/github/issues/ShanJiangEmugen/gene-target-screener" alt="Issues">
  <img src="https://img.shields.io/badge/code%20style-black-000000" alt="Code Style: black">
</p>


A simple and flexible tool for scanning, extracting, and aligning subsequences between two nucleotide sequences.      
It provides a sliding-window–based subsequence generator and a Biopython PairwiseAligner wrapper to evaluate conservation, similarity, and alignment scores across two input sequences.      

This tool is suitable for general sequence comparison tasks, conserved-region discovery, and lightweight pairwise alignment workflows.      

## Features

- **Sliding-window extraction** of subsequences from any input DNA/RNA sequence    
- **Automatic filtering** of low-complexity CTG repeat–rich regions    
- **Pairwise alignment** using Biopython `PairwiseAligner`   
- **Conservation-based filtering** across any two sequences    
- **Batch mode** for running multiple sequence comparisons from a CSV file    
- Outputs structured **CSV reports** for downstream analysis    
- Works as both a **Python module** and a **CLI tool**


## Installation

```bash
git clone https://github.com/<your-username>/u7-target-screener.git
cd u7-target-screener
pip install -r requirements.txt
```

## Input Format

The batch runner accepts a CSV file with the following columns:        
| Column         | Description                                  |
| -------------- | -------------------------------------------- |
| `Comparison #` | Each comparison must contain exactly 2 rows. |
| `Region`       | Name/label for each sequence.                |
| `Sequence`     | Raw nucleotide sequence (A/C/G/T).           |


### Example:
```csv
Comparison #,Region,Sequence
1,SeqA,ACGTACGTACGT...
1,SeqB,TACGATCGATCG...
2,Human_exon5,ACTG...
2,Mouse_exon5,ACTG...
```

## Usage
Batch Mode

Run all pairwise comparisons in a CSV:    
```bash
python -m gene_target_screener.batch_run \
    -i examples/screening_example.csv \
    -o results \
    -t 0.8
```

This generates a directory:    
```
results/
└── pair1_SeqA_SeqB/
    ├── input_SeqA_refseq_SeqB.csv
    └── input_SeqB_refseq_SeqA.csv
```

Each output CSV includes:
- `Position`     
- `Sequence`    
- `Length`    
- `Matched`    
- `Conservation`    
- `Score`    
- `Score Max`    
- `Score-wise Conservation`    


## Python API Example
```python
from gene_target_screener.aligner import SubsequenceAligner
from Bio.Seq import Seq

seq1 = Seq("ACGTACGTACGT...")
seq2 = Seq("TACGATCGATCG...")

aligner = SubsequenceAligner(seq1, seq2)
aligner.make_query(seq_low=18, seq_high=36)
aligner.init_aligner()
aligner.get_alignment(threshold=0.8)

df = aligner.results
print(df.head())
```

## Project Structure
```
gene-target-screener/
├── gene_target_screener/
│   ├── __init__.py
│   ├── aligner.py
│   └── batch_run.py
├── examples/
│   ├── screening_example.csv
│   └── run_example.sh
├── notebooks/
│   └── batch_run.ipynb
├── README.md
├── requirements.txt
└── LICENSE
```

## Requirements
```
biopython>=1.79
pandas>=1.4
numpy>=1.22
```

## License 
MIT      






