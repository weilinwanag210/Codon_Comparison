#!/usr/bin/env python
# Codon_Comparison.py
import os
import argparse
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
import csv

def codon_alignment_pretty_print_with_score_and_report(
    ref_id, ref_cds, var_id, var_cds, aligned_prot1, aligned_prot2, score_
):
    i1 = i2 = 0
    codon_line1, codon_line2, codon_line3 = [], [], []
    mutation_records = []
    total_score = 0

    for pos, (aa1, aa2) in enumerate(zip(aligned_prot1, aligned_prot2), start=1):
        codon1 = "---" if aa1 == '-' else ref_cds[i1*3:(i1+1)*3]; i1 += (aa1 != '-')
        codon2 = "---" if aa2 == '-' else var_cds[i2*3:(i2+1)*3]; i2 += (aa2 != '-')
        marker = "".join(['|' if c1 == c2 and aa1 == aa2 and c1 != "-" else '.' for c1, c2 in zip(codon1, codon2)])


        # scoring + classification
        if codon1 == codon2 and codon1 != "NNN":
            score = score_[0]
            mtype = "same"
        elif codon1 or codon2 == "NNN":
            score = score_[1]
            mtype = "NNN"
        elif codon1 == "---" or codon2 == "---":
            score = score_[2]
            mtype = "gap"
        elif aa1 == aa2 and aa1 != "-":
            score = score_[3]
            mtype = "synonymous"
        else:
            score = score_[4]
            mtype = "nonsynonymous"

        total_score += score

        # collect per-site mutation
        if mtype != "same":
            mutation_records.append({
                "Position": pos,
                "RefCodon": codon1,
                "VarCodon": codon2,
                "RefAA": aa1,
                "VarAA": aa2,
                "Type": mtype,
                "Score": score
            })

        # format for alignment view
        if aa1 == aa2 and aa1 != "-":
            codon_line1.append(f"({codon1} ")
            codon_line2.append(f" {marker} ")
            codon_line3.append(f" {codon2})")
        else:
            codon_line1.append(codon1)
            codon_line2.append(marker)
            codon_line3.append(codon2)

    # format alignment view block-wise
    block_size = 20
    id_width = 6
    output = []
    for i in range(0, len(codon_line1), block_size):
        blk1 = codon_line1[i:i+block_size]
        blk2 = codon_line2[i:i+block_size]
        blk3 = codon_line3[i:i+block_size]
        output.append(f"{ref_id:<{id_width}}: {' '.join(blk1)}")
        output.append(f"{'':<{id_width}}  {' '.join(blk2)}")
        output.append(f"{var_id:<{id_width}}: {' '.join(blk3)}")
        output.append("")

    output.append(f"{'Score':<{id_width}}: Codon alignment score = {total_score}")
    return "\n".join(output), mutation_records

def run_codon_alignment(cds_file, protein_file, output_file, report_file, score_, pre_aligned):
    with open(cds_file) as f1, open(protein_file) as f2:
        cds_iter = SeqIO.parse(f1, "fasta")
        prot_iter = SeqIO.parse(f2, "fasta")

        ref_cds_rec = next(cds_iter)
        ref_prot_rec = next(prot_iter)
        ref_id = ref_cds_rec.id
        ref_cds = str(ref_cds_rec.seq)
        ref_prot = str(ref_prot_rec.seq)

        with open(output_file, "w") as fout, open(report_file, "w", newline='') as ftable:
            writer = csv.DictWriter(ftable, fieldnames=["VarID", "Position", "RefCodon", "VarCodon", "RefAA", "VarAA", "Type", "Score"], delimiter='\t')
            writer.writeheader()
            if pre_aligned:
                for var_cds_rec, var_prot_rec in zip(cds_iter, prot_iter):
                    var_id = var_cds_rec.id
                    var_cds = str(var_cds_rec.seq)
                    var_prot = str(var_prot_rec.seq)
                    aln_result, mutation_table = codon_alignment_pretty_print_with_score_and_report(
                        ref_id, ref_cds, var_id, var_cds, ref_prot, var_prot, score_
                    )
                    fout.write(aln_result + "\n\n")
                    for row in mutation_table:
                        row["VarID"] = var_id
                        writer.writerow(row)
            else:        
                for var_cds_rec, var_prot_rec in zip(cds_iter, prot_iter):
                    var_id = var_cds_rec.id
                    var_cds = str(var_cds_rec.seq)
                    var_prot = str(var_prot_rec.seq)

                    align = pairwise2.align.globalms(ref_prot, var_prot, *score_[:4])[0]
                    aligned_prot1, aligned_prot2 = align.seqA, align.seqB

                    aln_result, mutation_table = codon_alignment_pretty_print_with_score_and_report(
                        ref_id, ref_cds, var_id, var_cds, aligned_prot1, aligned_prot2, score_
                    )
                    fout.write(aln_result + "\n\n")
                    for row in mutation_table:
                        row["VarID"] = var_id
                        writer.writerow(row)

    print(f" Codon alignment saved to {output_file}\n Mutation table saved to {report_file}")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Codon-level alignment with scoring and mutation reporting.")
    parser.add_argument("--cds", required=True, help="CDS aligned fasta file")
    parser.add_argument("--protein", required=True, help="Protein aligned fasta file")
    parser.add_argument("--output", required=True, help="Output alignment file")
    parser.add_argument("--report", required=True, help="TSV file to write mutation table")
    parser.add_argument("--score", nargs=5, type=int, default=[2, 0, -1, 1, -2],
                        metavar=("SAME", "NNN", "GAP", "SYN", "NONSYN"),
                        help="Five scores: same, NNN, gap, synonymous, nonsynonymous")
    parser.add_argument(
    "--pre_aligned", action="store_true",
    help="Use this flag if input CDS/protein files are already aligned"
)
    args = parser.parse_args()
    run_codon_alignment(args.cds, args.protein, args.output, args.report, args.score, args.pre_aligned)

