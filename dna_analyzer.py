#!/usr/bin/env python3
"""
DNA Analysis Toolkit
====================

A bioinformatics tool for DNA sequence analysis combining custom algorithms
with Biopython integration.

Author: Karolina Pajdak
Portfolio Project: July 2025
GitHub: github.com/karo-pajdak/DNA-sequence-analysis-tool

Core Features:
- DNA sequence validation and cleaning
- Reverse complement calculation
- GC content analysis
- DNA to protein translation
- ORF (Open Reading Frame) detection
- Restriction enzyme site mapping
- FASTA file processing
- Results export (JSON/CSV)

Usage: update this
    python dna_analyzer.py --sequence ATCGATCG --all
    python dna_analyzer.py --file sequences.fasta --find-orfs
    python dna_analyzer.py --help

Requirements: Python 3.7+, Biopython
"""
from Bio.Data import CodonTable

def reverse_complement(sequence):
    complement = []
    for base in sequence:
        if base == "A":
            complement.append("T")
        elif base == "T":
            complement.append("A")
        elif base == "C":
            complement.append("G")
        elif base == "G":
            complement.append("C")
    reverse_comp = ''.join(complement[::-1])
    return reverse_comp

def gc_content(sequence):
    count = 0
    for base in sequence:
        if base == "G" or base == "C":
            count = count + 1
    gc = count / len(sequence) * 100
    return gc

def translate(sequence):
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    start_codon_index = sequence.find("ATG")
    stop_codon = standard_table.stop_codons
    translated = ""
    if start_codon_index > -1:
        for nt in range(start_codon_index, len(sequence) - 2, 3):
            codon = sequence[nt:nt+3]
            if codon in stop_codon:
                break
            amino_acid = standard_table.forward_table.get(codon, "X")
            translated += amino_acid
    else:
        translated = "No start codon."
    return translated

def find_orfs(sequence):
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    start_codon = standard_table.start_codons
    stop_codon = standard_table.stop_codons
    ORFs = []
    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            if codon in start_codon:
                for j in range(i+3,len(sequence)-2,3):
                    stop = sequence[j:j+3]
                    if stop in stop_codon:
                        ORF_sequence = sequence[i:j+3]
                        ORFs.append(ORF_sequence)
                        break
                i += 3
            else:
                i += 3
    return ORFs


def analyze_sequence(sequence):
    print(f"Reverse Complement: {reverse_complement(sequence)}")
    print(f"GC Content: {gc_content(sequence)}")
    print(f"Translated sequence: {translate(sequence)}")
    print(f"ORFs: {find_orfs(sequence)}")
    

if __name__ == "__main__":
    test_dna = "AGCGTTGATGCAGTGCGTTGGTACGA"
    analyze_sequence(test_dna)