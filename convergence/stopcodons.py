#!/usr/bin/python3

# Goal: Convert aligned FASTA to PHYLIP for PAML.
# Robustly replace ANY real stop codon with 'NNN' while preserving gaps.

import os
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def robust_replace_stops_with_NNN(seq):
    """
    Replace real stop codons with NNN while preserving gaps.
    Handles ambiguous bases, non-standard characters, and
    sequences not divisible by 3.
    """
    STOP_CODONS = {'TAA', 'TAG', 'TGA'}

    seq = str(seq).upper()

    # Replace non-standard characters with N first
    clean_chars = set('ATCGNRYSWKMBDHV-')
    seq = ''.join(c if c in clean_chars else 'N' for c in seq)

    new_seq = list(seq)
    bases = []
    base_pos = []

    for i, char in enumerate(seq):
        if char != '-':
            bases.append(char)
            base_pos.append(i)

    # Only walk complete codons — ignore trailing partial codon
    n_complete_codons = len(bases) // 3

    for codon_idx in range(0, n_complete_codons * 3, 3):
        codon = ''.join(bases[codon_idx:codon_idx+3])

        # Only replace unambiguous stop codons
        # Skip codons with N or ambiguous bases — can't reliably call stop
        if codon in STOP_CODONS:
            for p in base_pos[codon_idx:codon_idx+3]:
                new_seq[p] = 'N'
        elif 'N' in codon or not all(c in 'ATCG' for c in codon):
            # Ambiguous codon — replace with NNN to be safe
            for p in base_pos[codon_idx:codon_idx+3]:
                new_seq[p] = 'N'

    return Seq(''.join(new_seq))

def fasta_to_phylip(aligned_fastas, output_phylip):
    # Make sure output directory exists
    if not os.path.exists(output_phylip):
        os.makedirs(output_phylip)

    # Gather input FASTA files
    input_files = [f for f in os.listdir(aligned_fastas) if f.endswith('_NT.fa') and not f.startswith(".")]

    for file in input_files:
        path = os.path.join(aligned_fastas, file)

        try:
            # Read alignment
            with open(path, "r", encoding='latin-1') as p:
                alignment = AlignIO.read(p, "fasta")

                # Process each sequence robustly
                new_alignment = MultipleSeqAlignment([
                    SeqRecord(
                        robust_replace_stops_with_NNN(record.seq),
                        id=record.id[:7],
                        description=""
                    )
                    for record in alignment
                ])

                # Optional sanity check
                for rec in new_alignment:
                    aa = Seq(str(rec.seq).replace("-", "")).translate(table=1)
                    if "*" in aa:
                        print(f"WARNING: {rec.id} still has stop codon after cleaning!")

                # Write PHYLIP
                output_file = os.path.join(output_phylip, f"{os.path.splitext(file)[0]}.phy")
                with open(output_file, "w") as out_handle:
                    AlignIO.write(new_alignment, out_handle, "phylip-sequential")

                print(f"Converted {file} -> PHYLIP (all stop codons robustly masked with NNN).")

        except Exception as e:
            print(f"Error processing {file}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert aligned FASTA files to PHYLIP format, replacing ANY real stop codon with NNN while preserving gaps."
    )
    parser.add_argument("aligned_fastas", help="Folder with aligned FASTA files")
    parser.add_argument("output_phylip", help="Output folder for PHYLIP files")
    args = parser.parse_args()

    fasta_to_phylip(args.aligned_fastas, args.output_phylip)
