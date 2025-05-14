#!/usr/bin/env python

import os
import random
import argparse
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment

CODON_TABLE = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "X": ["NNN"]
}

def simulate_variant(seq, x_rate=0.03, gap_rate=0.03, mutation_rate=0.02):
    AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
    result = []
    for aa in seq:
        r = random.random()
        if r < gap_rate:
            result.append('')
        elif r < gap_rate + x_rate:
            result.append(aa)
            result.extend(['X'] * random.choice([5, 12]))
        elif r < gap_rate + x_rate + mutation_rate:
            alt = random.choice([a for a in AA_ALPHABET if a != aa])
            result.append(alt)
        else:
            result.append(aa)
   
    cds = "".join(random.choice(CODON_TABLE[aa]) for aa in result if aa in CODON_TABLE)
    return ''.join(result), cds

def generate_simulated_alignment(ref_seq, num_variants=5, output_dir="simulated_output"):
    os.makedirs(output_dir, exist_ok=True)
    ref_seq = ref_seq.replace("\n", "")

    prot_fasta = os.path.join(output_dir, "protein_aligned.fasta")
    cds_fasta = os.path.join(output_dir, "cds_aligned.fasta")
    txt_output = os.path.join(output_dir, "alignment_view.txt")

    ref_prot = SeqRecord(Seq(ref_seq), id="Ref", description="Reference")
    ref_cds = SeqRecord(Seq("".join(random.choice(CODON_TABLE[a]) for a in ref_seq)), id="Ref", description="CDS")
    prot_records = [ref_prot]
    cds_records = [ref_cds]

    with open(txt_output, "w") as fout:
        for i in range(num_variants):
            prot, cds = simulate_variant(ref_seq)
            aln = pairwise2.align.globalms(ref_seq, prot, 2, -1, -10, -0.5)[0]
            aligned_prot = aln.seqB
            prot_records.append(SeqRecord(Seq(aligned_prot), id=f"Var{i+1}", description=''))
            cds_records.append(SeqRecord(Seq(cds), id=f"Var{i+1}", description=''))
            fout.write(f">>> Variant {i+1}\n")
            fout.write(format_alignment(*aln) + "\n")

    SeqIO.write(prot_records, prot_fasta, "fasta")
    SeqIO.write(cds_records, cds_fasta, "fasta")
    print(f"\n✅ Protein alignment saved to: {prot_fasta}")
    print(f"✅ CDS sequences saved to:     {cds_fasta}")
    print(f"📄 Alignment log saved to:     {txt_output}")

def main():
    parser = argparse.ArgumentParser(description="Simulate protein variants and CDSs from a given FASTA reference.")
    parser.add_argument("--ref", type=str, required=True, help="Input FASTA file of reference protein sequences")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--n", type=int, default=5, help="Number of variants to simulate per reference")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (optional)")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    os.makedirs(args.outdir, exist_ok=True)
    for i, record in enumerate(SeqIO.parse(args.ref, "fasta"), start=1):
        ref_name = record.id
        ref_seq = str(record.seq).replace("\n", "")
        subdir = os.path.join(args.outdir, f"{ref_name}")
        print(f"\n▶ Generating variants for {ref_name} into {subdir}")
        generate_simulated_alignment(ref_seq, num_variants=args.n, output_dir=subdir)

if __name__ == "__main__":
    main()
