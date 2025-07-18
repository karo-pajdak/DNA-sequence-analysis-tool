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

Usage: 
    python3 dna_analyzer.py --sequence ATCGATCG --all
    python3 dna_analyzer.py --file sequences.fasta --find-orfs
    python3 dna_analyzer.py --file sequences.fasta --restriction_enzymes EcoRI HaeIII
    python3 dna_analyzer.py --help

Requirements: Python 3.7+, Biopython
"""
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, CommOnly
import os
import argparse

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

def read_file(filename):
    extension = os.path.splitext(filename)[1].lower()
    sequences = []
    if extension in [".fasta", ".fa", ".fna", ".ffn"]:
        format_type = "fasta"
    elif extension in [".fastq", ".fq"]:
        format_type = "fastq"
    else:
        raise ValueError("Unsupported file format.")
    for record in list(SeqIO.parse(filename, format_type)):
        sequences.append(str(record.seq))
    if not sequences:
        print("No sequences found in file.")
        return []
    return sequences

def find_restriction_sites(sequence, restriction_enzymes):
    valid_enzymes = {enzyme.__name__ for enzyme in CommOnly}
    invalid_user_enzymes = [enz for enz in restriction_enzymes if enz not in valid_enzymes]
    valid_user_enzymes = [enz for enz in restriction_enzymes if enz in valid_enzymes]
    if invalid_user_enzymes:
        print("Invalid enzymes below will be skipped:")
        for enz in invalid_user_enzymes:
            print(f"{enz}")
    if not valid_user_enzymes:
        print("No valid enzymes provided.")
        return {}
    enzyme_classes = [e for e in CommOnly if e.__name__ in valid_user_enzymes]
    restriction_batch = RestrictionBatch(enzyme_classes)
    restriction_sites = restriction_batch.search(sequence)
    return restriction_sites

def valid_sequence(sequence):
    valid_bases = set("ATGCN")
    sequence = sequence.upper()
    return all(base in valid_bases for base in sequence)

def main():
    parser = argparse.ArgumentParser(description="The DNA sequence analysis tool allows users to find the reverse complement, GC content, sequence translation, ORF, and restriction enzymes from either a command line sequence or FASTA or FASTQ file.")

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--file', help='Path to input FASTA or FASTQ file.')
    input_group.add_argument('--sequence', help='DNA sequence string.')

    parser.add_argument('--reverse_complement', action='store_true', help='Get reverse complement')
    parser.add_argument('--gc_content', action='store_true', help='Get GC content')
    parser.add_argument('--translate', action='store_true', help='Translate DNA to protein')
    parser.add_argument('--find_orfs', action='store_true', help="Find open reading frames")
    parser.add_argument('--restriction_enzymes', nargs='+', metavar='ENZYME', help="Find cut sites for given restriction enzymes (ex: EcoRI)")
    parser.add_argument('--all', action='store_true', help="Perform all actions")

    args = parser.parse_args()

    if args.all:
        args.reverse_complement = True
        args.gc_content = True
        args.translate = True
        args.find_orfs = True
        if not args.restriction_enzymes:
            from Bio.Restriction import CommOnly
            args.restriction_enzymes = [e.__name__ for e in CommOnly]

    if args.sequence:
        sequences = [args.sequence]
    else:
        sequences = read_file(args.file)

    for i, seq in enumerate(sequences):
        sequence = seq.upper().replace('\n', '').replace(' ', '')
        print(f"Sequence {i+1}: {sequence}")
        if not valid_sequence(sequence):
            print(f"Sequence {i} contains invalid characters. Skipped.")
            continue
        if args.reverse_complement:
            print("Reverse complement: ", reverse_complement(sequence))
        if args.gc_content:
            print("GC content: {:.2f}%".format(gc_content(sequence)))
        if args.translate:
            print("Translated sequence: ", translate(sequence))
        if args.find_orfs:
            orfs = find_orfs(sequence)
            print(f"Found {len(orfs)} ORFs")
            for j, orf in enumerate(orfs, 1):
                print(f"ORF {j}: {orf}")
        if args.restriction_enzymes:
            enzyme_results = find_restriction_sites(Seq(sequence), args.restriction_enzymes)
            sites_found = {enzyme: pos for enzyme, pos in enzyme_results.items() if pos}
            if sites_found:
                print(f"Restriction enzyme cut sites found:")
                for enzyme, positions in sites_found.items():
                    print(f"{enzyme}: {positions}")
            else:
                print("No restriction enzyme sites found.")
        print()


if __name__ == "__main__":
    main()