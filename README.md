# DNA-sequence-analysis-tool
A bioinformatics tool for DNA sequence analysis combining custom algorithms with Biopython integration.

This tool performs key DNA analyses like reverse complementing, GC content calculation, translation to protein, ORF detection, and restriction enzyme site mapping. It supports both direct sequence input and file input (`.fasta` or `.fastq`).

## Features
- DNA sequence validation and cleaning
- Reverse complement calculation
- GC content analysis
- DNA to protein translation
- ORF (Open Reading Frame) detection
- Restriction enzyme site mapping
- FASTA file processing

## Folder Structure
```
DNA-sequence-analysis-tool/
├─ dna_analyzer.py # Main Python script for analysis
├─ README.md # Project documentation and instructions
└─ .gitignore # Git ignore file to exclude large files (ex: FASTA, FASTQ)
```

## Quick Start
1. Install Biopython
```
pip install biopython
```

2. Run the scrript
```
python3 dna_analyzer.py --sequence ATGCGTACGTAGCT --all
```

## Example Commands
## Analyze a single DNA sequence
```
python3 dna_analyzer.py --sequence ATGCGTACGTAGCT --reverse_complement --gc_content
```
### Process sequences from a FASTA file
```
python3 dna_analyzer.py --file my_sequence.fasta --translate
```
### Find open reading frames and restriction enzyme sites for a specific enzyme
```
python3 dna_analyzer.py --file my_sequence.fasta --find_orfs --restriction_enzymes EcoRI
```
### Run all available analyses
```
python3 dna_analyzer.py --file my_sequences.fasta --all
```

## Supported File Types
- `.fasta,` `.fa`, `.fna`, `.ffn`
- `.fastq`, `.fq`


## Arguments Summary
| Flag                  | Description                                               |
|-----------------------|-----------------------------------------------------------|
| `--file`              | Input FASTA/FASTQ file                                    |
| `--sequence`          | Direct DNA sequence input                                 |
| `--reverse_complement`| Get the reverse complement                               |
| `--gc_content`        | Calculate GC%                                            |
| `--translate`         | Translate to protein                                     |
| `--find_orfs`         | Find open reading frames                                 |
| `--restriction_enzymes` | Search for restriction enzyme sites (e.g. EcoRI, BamHI) |
| `--all`               | Run all features on all sequences                        |


## Author
### Karolina Pajdak
Bioinformatics Portfolio Project – July 2025

github.com/karo-pajdak
