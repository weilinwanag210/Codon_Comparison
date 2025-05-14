# Codon_Comparison
DEMO For Codon Comparison
## Quick Start
Setting Python Environment and dependent packages directly

`conda env create -f environment.yml`
## Test Data Generation
### **Generated Variant Sequences based on Given Reference Protein Data**
```python
ref_seq = (
    "MMARPVDPQRSPDPTFRSSTRHSGKLEPMEATAHLLRKQCPSRLNSPAWEASGLHWSSLDSPVGSMQALRPSAQHSWS"
    "PEPSVVPDQAWEDTALHQKKLCPLSLTSLPREAAVNFSYRSQTLLQEAQVLQGSPELLPRSPKPSGLQRLAPEEATAL"
    "PLRRLCHLSLMEKDLGTTAHPRGFPELSHKSTAAASSRQSRPRVRSASLPPRTRLPSGSQAPSAAHPKRLSDLLLTSRA"
    "AAPGWRSPDPRSRLAAPPLGSTTLPSTWTAPQSRLTARPSRSPEPQIRESEQRDPQLRRKQQRWKEPLMPRREEKYPLR"
    "GTDPLPPGQPQRIPLPGQPLQPQPILTPGQPQKIPTPGQHQPILTPGHSQPIPTPGQPLPPQPIPTPGRPLTPQPIPTP"
    "GRPLTPQPIQMPGRPLRLPPPLRLLRPGQPMSPQLRQTQGLPLPQPLLPPGQPKSAGRPLQPLPPGPDARSISDPPAPR"
    "SRLPIRLLRGLLARLPGGASPRAAAAAACTTMKGWPAATMTPAETSPTMGPPDASAGFSIGEIAAAESPSATYSATFSC"
    "KPSGAASVDLRVPSPKPRALSRSRRYPWRRSADRCAKKPWRSGPRSAQRRNAVSSSTNNSRTKRWATCVRTACCF"
)
```
### **There are 3 types of mutation: Gap, Replace, Vague Insertion**
One aa in reference will miss or be replaced by another amino acid or be inserted *XXXXX*in variant sequence with a fixed probablity
```python
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
```
### Generated CDS based on the reference and generated variant aa sequences
```python
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
ref_cds = SeqRecord(Seq("".join(random.choice(CODON_TABLE[a]) for a in ref_seq)), id="Ref", description="CDS")
CDS = "".join(
    random.choice(CODON_TABLE[aa])
    for aa in result
    if aa in CODON_TABLE
    )
```
### Multiple Sequence Alignment

Call Mafft to align reference and variant sequences

`bash inputs.sh &&     mafft --auto          --thread ${GALAXY_SLOTS:-1} --threadit 0  input.fa > '/jetstream2/scratch/main/jobs/67636675/outputs/dataset_a6b3f507-e4aa-47e4-8ad7-73754666f940.dat'`
 
## Codon Comparison
Input cds.fasta & mafft_aligned_protein.fasta
`python Codon_Comparison.py \
  --cds cds_aligned.fasta \
  --protein protein_aligned.fasta \
  --output codon_alignment_results.txt \
  --report mutation_table.tsv \
  --score 2 0 -1 1 -2 \
  --pre_aligned
`
if input aa sequences have been aligned in advance please add `--pre_alligned`
### Output
The output will contain alignment file and mutantation file.
Alignment file looks as follows:
    Codons enclosed in parentheses indicate synonymous amino acids. The presence of a synonymous mutation is marked with `|` if the codons differ but encode the same amino acid.
    Scoring scheme for codon alignment:
    Match: +2
    Synonymous mutation: +1
    Gap: –1
    Nonsynonymous mutation: –2
    NNN: 0
```python
ref:  (ATG  ATG (GCC  (CGA  (CCT  (GTT  (GAT  (CCG  (CAG  (CGC  (TCA  (CCT  (GAT  (CCT  (ACG  (TTC  CGA TCA (AGT  (ACA  (CGG  (CAC  (AGC  GGG AAG (TTA  (GAG  (CCA  (ATG  (GAG  (GCT  (ACA  (GCT  (CAC  (TTA  --- --- --- --- --- (CTT  CGT (AAA  (CAA  TGC (CCA  (TCT  (AGG  (CTT  (AAT  (TCG  (CCA  (GCG  (TGG  (GAA  (GCG  (AGC  (GGT  (CTA  (CAT  (TGG  (TCA  (AGT  (CTC  (GAC  (AGT  (CCG  (GTA  (GGT  (AGT  (ATG  (CAA  (GCG  (TTG  --- --- --- --- --- (CGT  (CCC  (TCA  (GCA  (CAA  (CAT  (TCA  (TGG  (AGC  (CCT  GAA (CCT  (TCG  (GTG  (GTC  CCA (GAT  (CAG  (GCG  (TGG  (GAG  (GAC  (ACC  (GCG  (CTG  (CAC  (CAG  (AAG  (AAA  (CTC  (TGT  (CCA  (CTT  (TCG  (TTG  (ACT  (AGC  (CTG  (CCG  (AGG  (GAA  (GCG  (GCA  (GTA  AAT (TTT  (AGC  (TAT  (CGG  (AGT  (CAA  (ACT  (CTG  (CTC  (CAA  (GAG  GCA (CAA  (GTC  (TTG  (CAG  (GGT  (TCG  (CCC  (GAG  (CTG  (TTA  (CCT  (CGT  (AGT  (CCG  (AAA  (CCC  (TCT  GGA (CTC  (CAG  --- --- --- --- --- (CGG  (CTA  (GCG  (CCG  (GAG  (GAA  (GCT  (ACT  (GCA  (CTC  CCA (TTA  --- --- --- --- --- --- --- --- --- --- --- --- (CGC  (AGG  (CTG  (TGT  (CAC  (CTA  (AGT  (CTT  (ATG  (GAA  (AAA  (GAC  (CTT  (GGC  (ACA  (ACG  (GCA  (CAC  (CCC  (AGG  (GGC  (TTC  (CCA  (GAA  (CTG  (TCC  (CAT  (AAG  (TCA  (ACG  (GCG  (GCA  (GCT  (TCT  (TCG  (CGT  (CAA  --- --- --- --- --- (TCG  (AGA  (CCG  (CGC  --- --- --- --- --- (GTA  (AGG  (TCA  (GCA  (TCT  --- --- --- --- --- --- --- --- --- --- CTG (CCT  (CCG  (CGC  (ACT  (CGA  (TTA  (CCG  (TCC  (GGC  (AGC  (CAG  (GCT  (CCG  (AGC  (GCT  (GCT  (CAT  (CCT  (AAA  (CGA  (CTA  (AGC  (GAC  (CTC  TTG (TTG  (ACT  (TCA  (CGC  (GCC  (GCC  (GCT  (CCA  (GGA  (TGG  (CGA  (AGT  (CCG  (GAT  (CCG  (CGT  (TCT  (AGG  (CTG  (GCA  (GCT  (CCG  (CCG  (CTG  (GGA  (TCC  (ACC  (ACG  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- TTG CCT (TCA  (ACC  (TGG  --- --- --- --- --- (ACT  (GCC  (CCA  (CAA  (TCT  (CGG  (TTA  (ACA  (GCG  (CGA  (CCG  (TCC  (CGA  (TCG  (CCC  (GAA  (CCA  (CAA  (ATC  (CGA  (GAG  (AGC  --- --- --- --- --- (GAG  (CAG  (AGG  (GAC  (CCA  (CAG  (CTA  (CGC  (CGC  (AAA  (CAA  (CAA  (CGA  (TGG  (AAG  GAA (CCT  (CTA  (ATG  (CCC  (AGG  (AGG  (GAG  (GAG  (AAA  (TAC  (CCG  CTC (CGG  (GGA  --- --- --- --- --- (ACT  (GAC  (CCT  (CTA  --- --- --- --- --- --- --- --- --- --- --- --- (CCG  (CCG  (GGC  (CAG  (CCG  (CAG  --- --- --- --- --- (CGC  (ATA  (CCT  (TTG  (CCC  (GGG  (CAA  (CCG  (CTC  (CAG  CCC (CAA  (CCA  (ATA  (CTT  (ACA  (CCA  (GGG  (CAG  (CCG  (CAG  (AAG  (ATA  (CCG  (ACA  (CCT  (GGA  (CAA  (CAC  (CAA  (CCC  (ATT  --- --- --- --- --- (CTG  (ACA  (CCG  (GGA  (CAT  (TCA  (CAA  (CCG  (ATT  (CCG  (ACC  (CCT  (GGT  CAA (CCC  CTC (CCA  CCC (CAG  (CCG  (ATA  (CCG  (ACA  (CCC  (GGT  (AGA  (CCC  (CTT  (ACT  (CCC  (CAG  (CCA  (ATA  (CCG  (ACC  --- --- --- --- --- --- --- --- --- --- --- --- (CCG  (GGA  (CGA  (CCG  (CTG  (ACT  (CCG  (CAG  (CCG  (ATA  (CAG  (ATG  (CCC  (GGT  (AGG  (CCA  (TTG  (AGG  (TTA  (CCC  (CCG  (CCT  (TTG  (CGG  (CTA  (CTC  (CGC  (CCG  (GGC  (CAG  (CCC  (ATG  (AGC  (CCC  (CAG  (TTG  (CGT  (CAA  (ACG  (CAG  (GGC  (CTC  (CCC  TTG (CCT  (CAA  (CCC  (CTA  (CTG  (CCG  (CCA  (GGG  (CAG  (CCC  (AAG  (TCC  (GCT  GGT (AGG  (CCC  (CTC  (CAA  (CCG  (CTG  (CCG  (CCA  (GGT  (CCG  (GAC  (GCT  (CGT  (TCC  (ATA  (TCC  (GAC  (CCC  (CCG  (GCC  (CCA  (CGG  (AGC  (CGA  (TTA  (CCC  (ATC  (CGA  (CTG  (TTG  (AGA  (GGG  (CTT  CTG GCG (CGG  (CTA  (CCG  (GGC  (GGA  (GCC  (AGT  (CCG  (CGA  (GCG  (GCC  (GCC  (GCT  (GCG  (GCT  (TGT  (ACC  (ACA  (ATG  (AAG  --- --- --- --- --- (GGA  (TGG  (CCG  (GCC  (GCT  (ACC  (ATG  ACA (CCT  (GCT  --- --- --- --- --- --- --- --- --- --- --- --- (GAG  (ACG  (AGC  CCA (ACC  (ATG  (GGA  (CCT  (CCA  (GAC  (GCC  (TCA  (GCT  (GGC  (TTT  (AGC  (ATT  --- --- --- --- --- (GGC  (GAG  (ATC  (GCG  (GCA  (GCG  GAA (TCG  (CCT  (TCA  (GCC  (ACA  (TAT  (TCT  (GCT  (ACG  (TTT  (TCG  (TGT  (AAA  (CCA  (TCT  (GGA  (GCT  (GCC  (TCT  (GTA  (GAT  (TTG  (AGA  (GTT  (CCT  (TCG  (CCC  (AAA  (CCA  (AGA  (GCG  (CTT  (TCA  (AGG  (AGT  (AGG  (CGT  (TAT  (CCA  (TGG  (CGA  (CGT  (TCT  (GCT  (GAT  (AGG  (TGT  (GCT  --- --- --- --- --- --- --- --- --- --- --- --- (AAG  (AAA  (CCG  (TGG  (CGA  (TCA  (GGA  (CCC  (CGC  (TCT  (GCG  (CAG  (CGT  CGT (AAT  (GCC  (GTC  (AGT  (TCA  (AGT  (ACC  (AAC  (AAT  (AGT  CGG (ACC  (AAG  (CGT  (TGG  (GCC  (ACG  (TGT  (GTC  (CGT  (ACT  (GCT  (TGC  (TGT  (TTC 
       |||  ...  |||   ||.   ||.   |||   ||.   ||.   |||   ||.   ||.   |||   ||.   ||.   ||.   ||.  ... ...  ||.   |||   ||.   ||.   ...  ... ...  ||.   ||.   ||.   |||   ||.   ||.   ||.   ||.   |||   .|.  ... ... ... ... ...  ||.  ...  |||   |||  ...  ||.   ...   .|.   .|.   ||.   ...   ||.   |||   |||   |||   ||.   ...   ||.   |||   |||   |||   ||.   |||   .|.   ||.   ...   ||.   |||   ||.   ||.   |||   ||.   |||   .||  ... ... ... ... ...  ||.   ||.   ||.   ||.   |||   ||.   ...   |||   ...   |||  ...  ||.   ||.   |||   |||  ...  ||.   |||   ||.   |||   ||.   ||.   ||.   ||.   ||.   ||.   ||.   |||   |||   .|.   ||.   ||.   ||.   ||.   .|.   |||   ||.   ||.   ||.   .|.   ||.   ||.   ||.   ||.  ...  |||   ..|   ||.   .|.   ...   ||.   |||   ||.   .|.   |||   ||.  ...  |||   ||.   ||.   ||.   |||   |||   ||.   |||   ||.   .|.   ||.   ||.   ..|   ||.   ||.   ||.   ..|  ...  ||.   |||  ... ... ... ... ...  ||.   |||   ||.   ||.   |||   |||   ||.   ||.   ||.   |||  ...  .|.  ... ... ... ... ... ... ... ... ... ... ... ...  .|.   .|.   ||.   ||.   ||.   ||.   ...   |||   |||   |||   |||   ||.   |||   |||   ||.   |||   ||.   ||.   ||.   .|.   ||.   ||.   ||.   ||.   ||.   ...   |||   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   |||  ... ... ... ... ...  ||.   ||.   ||.   .|.  ... ... ... ... ...  ||.   .|.   ...   ||.   ...  ... ... ... ... ... ... ... ... ... ... ...  ||.   ||.   |||   ||.   .||   .||   ||.   |||   |||   ...   |||   ||.   ||.   ||.   |||   ||.   |||   |||   ||.   ||.   .|.   |||   ||.   |||  ...  ||.   ||.   ||.   .|.   ||.   ||.   ||.   ||.   |||   |||   .||   ...   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   |||   ||.   |||   ||.   ||.   ||.   ||.  ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ...  ...   ||.   |||  ... ... ... ... ...  ||.   ||.   |||   ||.   ||.   |||   .||   ||.   ||.   ||.   ||.   |||   ||.   ...   ||.   |||   ||.   |||   ||.   .|.   |||   ||.  ... ... ... ... ...  ||.   |||   |||   |||   |||   ||.   .|.   |||   ||.   ||.   ||.   ||.   ||.   |||   ||.  ...  ||.   ||.   |||   |||   .|.   ||.   |||   |||   |||   ||.   ||.  ...  ||.   |||  ... ... ... ... ...  ||.   |||   |||   .||  ... ... ... ... ... ... ... ... ... ... ... ...  ||.   ||.   |||   ||.   ||.   ||.  ... ... ... ... ...  ||.   |||   ||.   .|.   ||.   ||.   ||.   ||.   .|.   ||.  ...  |||   |||   ||.   ||.   ||.   ||.   |||   ||.   ||.   |||   ||.   ||.   ||.   ||.   ||.   ||.   ||.   |||   ||.   ||.   |||  ... ... ... ... ...  |||   ||.   ||.   |||   |||   |||   |||   ||.   ||.   ||.   ||.   |||   ||.  ...  ||.  ...  ||.  ...  ||.   |||   ||.   ||.   ||.   ||.   ||.   .||   |||   |||   ||.   |||   |||   |||   |||   |||   ||.  ... ... ... ... ... ... ... ... ... ... ... ...  ||.   |||   |||   ||.   ||.   ||.   ||.   |||   ||.   |||   |||   |||   |||   |||   .|.   ||.   .|.   |||   |||   ||.   ||.   ||.   .||   ||.   ||.   ||.   ||.   |||   ||.   |||   ||.   |||   ...   |||   ||.   .|.   |||   |||   |||   |||   ||.   ||.   ||.  ...  ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   |||   |||   ||.  ...  .||   |||   .|.   |||   ||.   ||.   |||   ||.   ||.   ||.   ||.   ||.   ||.   ..|   ||.   ...   ||.   ||.   ||.   ||.   ||.   ||.   |||   |||   .|.   ||.   ||.   ||.   ||.   .||   .||   |||   ||.  ... ...  .|.   ||.   ||.   ||.   ||.   |||   ..|   ||.   ||.   ||.   ||.   |||   ||.   |||   ||.   ||.   ||.   |||   |||   ||.  ... ... ... ... ...  ||.   |||   ||.   ||.   ||.   ||.   |||  ...  ||.   |||  ... ... ... ... ... ... ... ... ... ... ... ...  ||.   ||.   ||.  ...  ||.   |||   ||.   |||   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ...   |||  ... ... ... ... ...  ||.   ||.   |||   ||.   |||   ||.  ...  ...   ||.   ...   ||.   |||   |||   ||.   |||   ||.   |||   ...   ||.   |||   |||   ||.   |||   ||.   ||.   ||.   |||   |||   ||.   ||.   ||.   ||.   ||.   ||.   |||   ||.   .|.   ||.   ||.   ||.   ||.   ...   |||   |||   ||.   ||.   |||   .||   .|.   ..|   ||.   |||   ||.   |||   ||.  ... ... ... ... ... ... ... ... ... ... ... ...  |||   ||.   ||.   |||   .||   ||.   |||   ||.   ||.   |||   ||.   |||   |||  ...  ||.   ||.   ||.   ||.   ...   ...   ||.   ||.   |||   ||.  ...  ||.   ||.   .|.   |||   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||.   ||. 
var:   ATG) CAC  GCC)  CGT)  CCG)  GTT)  GAC)  CCT)  CAG)  CGA)  TCG)  CCT)  GAC)  CCC)  ACA)  TTT) CAG ---  AGC)  ACA)  CGT)  CAT)  TCG) GCA ---  TTG)  GAA)  CCG)  ATG)  GAA)  GCG)  ACG)  GCC)  CAC)  CTG) NNN NNN NNN NNN NNN  CTG) ---  AAA)  CAA) ACA  CCG)  AGC)  CGC)  TTA)  AAC)  AGC)  CCG)  GCG)  TGG)  GAA)  GCT)  TCT)  GGG)  CTA)  CAT)  TGG)  TCC)  AGT)  TTG)  GAT)  TCG)  CCT)  GTA)  GGG)  AGC)  ATG)  CAG)  GCG)  CTG) NNN NNN NNN NNN NNN  CGC)  CCG)  TCC)  GCC)  CAA)  CAC)  AGC)  TGG)  TCA)  CCT) ---  CCC)  TCA)  GTG)  GTC) GTC  GAC)  CAG)  GCA)  TGG)  GAA)  GAT)  ACA)  GCT)  CTC)  CAT)  CAA)  AAG)  AAA)  TTG)  TGC)  CCT)  CTC)  TCT)  CTT)  ACT)  AGT)  CTA)  CCT)  CGA)  GAG)  GCT)  GCC)  GTC) ---  TTT)  TCC)  TAC)  AGA)  TCG)  CAG)  ACT)  CTA)  TTG)  CAA)  GAA) TTG  CAA)  GTA)  TTA)  CAA)  GGT)  TCG)  CCA)  GAG)  CTT)  CTC)  CCA)  CGG)  TCT)  CCA)  AAG)  CCA)  AGT) ---  CTA)  CAG) NNN NNN NNN NNN NNN  CGC)  CTA)  GCT)  CCT)  GAG)  GAA)  GCG)  ACG)  GCC)  CTC) ---  CTC) NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN  AGA)  CGT)  CTA)  TGC)  CAT)  CTC)  TCC)  CTT)  ATG)  GAA)  AAA)  GAT)  CTT)  GGC)  ACG)  ACG)  GCT)  CAT)  CCA)  CGA)  GGA)  TTT)  CCC)  GAG)  CTC)  AGT)  CAT)  AAA)  TCG)  ACA)  GCA)  GCG)  GCA)  TCG)  TCC)  CGC)  CAA) NNN NNN NNN NNN NNN  TCA)  AGG)  CCA)  AGG) NNN NNN NNN NNN NNN  GTT)  CGA)  AGT)  GCT)  AGC) NNN NNN NNN NNN NNN CTT NNN NNN NNN NNN NNN  CCG)  CCT)  CGC)  ACA)  AGA)  CTA)  CCT)  TCC)  GGC)  TCG)  CAG)  GCG)  CCT)  AGT)  GCT)  GCG)  CAT)  CCT)  AAG)  CGG)  TTG)  AGC)  GAT)  CTC) TAT  TTA)  ACC)  TCG)  AGG)  GCT)  GCG)  GCG)  CCT)  GGA)  TGG)  AGA)  TCA)  CCT)  GAC)  CCA)  CGA)  TCC)  AGA)  CTT)  GCC)  GCG)  CCG)  CCA)  CTG)  GGC)  TCT)  ACT)  ACC) NNN NNN NNN NNN NNN --- CCG NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN  AGC)  ACA)  TGG) NNN NNN NNN NNN NNN  ACG)  GCA)  CCA)  CAG)  TCA)  CGG)  CTA)  ACC)  GCC)  CGC)  CCC)  TCC)  CGG)  AGT)  CCA)  GAA)  CCG)  CAA)  ATT)  AGG)  GAG)  AGT) NNN NNN NNN NNN NNN  GAA)  CAG)  AGG)  GAC)  CCA)  CAA)  TTG)  CGC)  CGT)  AAG)  CAG)  CAG)  CGT)  TGG)  AAA) ---  CCG)  CTT)  ATG)  CCC)  CGA)  AGA)  GAG)  GAG)  AAA)  TAT)  CCC) TCT  CGA)  GGA) NNN NNN NNN NNN NNN  ACA)  GAC)  CCT)  TTA) NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN  CCC)  CCC)  GGC)  CAA)  CCA)  CAA) NNN NNN NNN NNN NNN  CGT)  ATA)  CCA)  CTT)  CCG)  GGC)  CAG)  CCC)  TTA)  CAA) TGG  CAA)  CCA)  ATT)  CTG)  ACG)  CCC)  GGG)  CAA)  CCA)  CAG)  AAA)  ATT)  CCC)  ACT)  CCG)  GGT)  CAG)  CAC)  CAG)  CCG)  ATT) NNN NNN NNN NNN NNN  CTG)  ACT)  CCT)  GGA)  CAT)  TCA)  CAA)  CCC)  ATC)  CCT)  ACA)  CCT)  GGG) ---  CCT) ---  CCC) ACC  CAA)  CCG)  ATT)  CCC)  ACG)  CCG)  GGG)  CGA)  CCC)  CTT)  ACA)  CCC)  CAG)  CCA)  ATA)  CCG)  ACA) NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN  CCT)  GGA)  CGA)  CCT)  CTT)  ACA)  CCA)  CAG)  CCA)  ATA)  CAG)  ATG)  CCC)  GGT)  CGT)  CCT)  CTC)  AGG)  TTA)  CCT)  CCC)  CCA)  CTG)  CGA)  CTT)  CTG)  CGT)  CCG)  GGA)  CAG)  CCT)  ATG)  TCA)  CCC)  CAA)  CTT)  CGT)  CAA)  ACG)  CAG)  GGA)  CTG)  CCG) ---  CCG)  CAG)  CCG)  CTG)  CTA)  CCT)  CCT)  GGA)  CAA)  CCG)  AAG)  TCC)  GCA) ATT  CGG)  CCC)  TTG)  CAA)  CCA)  CTC)  CCG)  CCG)  GGG)  CCT)  GAT)  GCA)  CGA)  AGC)  ATC)  AGT)  GAT)  CCT)  CCA)  GCA)  CCG)  CGA)  AGC)  CGA)  CTT)  CCA)  ATA)  CGT)  CTA)  CTG)  CGA)  GGG)  CTG) GCT ---  AGA)  CTT)  CCA)  GGA)  GGT)  GCC)  TCT)  CCT)  CGG)  GCT)  GCG)  GCC)  GCC)  GCG)  GCG)  TGC)  ACA)  ACA)  ATG)  AAA) NNN NNN NNN NNN NNN  GGT)  TGG)  CCA)  GCT)  GCA)  ACG)  ATG) GTC  CCC)  GCT) NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN  GAA)  ACA)  AGT) ---  ACA)  ATG)  GGG)  CCT)  CCT)  GAT)  GCT)  TCG)  GCG)  GGG)  TTC)  TCG)  ATT) NNN NNN NNN NNN NNN  GGA)  GAA)  ATC)  GCC)  GCA)  GCT) GTT  AGT)  CCG)  AGT)  GCA)  ACA)  TAT)  TCC)  GCT)  ACT)  TTT)  AGT)  TGC)  AAA)  CCA)  TCA)  GGA)  GCA)  GCA)  TCC)  GTA)  GAT)  TTA)  AGG)  GTG)  CCC)  TCT)  CCG)  AAA)  CCT)  CGT)  GCC)  CTC)  TCT)  AGA)  TCC)  AGG)  CGT)  TAC)  CCT)  TGG)  AGA)  AGG)  AGT)  GCG)  GAT)  AGA)  TGT)  GCA) NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN NNN  AAG)  AAG)  CCT)  TGG)  AGA)  TCG)  GGA)  CCT)  CGG)  TCT)  GCT)  CAG)  CGT) ---  AAC)  GCT)  GTG)  AGC)  AGC)  TCA)  ACG)  AAT)  AAT)  AGC) AAC  ACT)  AAA)  AGA)  TGG)  GCT)  ACC)  TGC)  GTG)  CGA)  ACA)  GCC)  TGT)  TGC)  TTT)

 Codon alignment score: 569
```
### Mutantion File output example
| VarID | Position | RefCodon | VarCodon | RefAA | VarAA | Type           | Score |
|-------|----------|----------|----------|--------|--------|----------------|--------|
| Var1  | 2        | ATG      | CAC      | M      | H      | nonsynonymous  | -2     |
| Var1  | 4        | CGA      | CGT      | R      | R      | synonymous     | 1      |
| Var1  | 5        | CCT      | CCG      | P      | P      | synonymous     | 1      |
| Var1  | 7        | GAT      | GAC      | D      | D      | synonymous     | 1      |
| Var1  | 8        | CCG      | CCT      | P      | P      | synonymous     | 1      |
| Var1  | 15       | CGC      | CGA      | R      | R      | synonymous     | 1      |
| Var1  | 16       | TCA      | TCG      | S      | S      | synonymous     | 1      |
| Var1  | 18       | GAT      | GAC      | D      | D      | synonymous     | 1      |
| Var1  | 19       | CCT      | CCC      | P      | P      | synonymous     | 1      |
| Var1  | 20       | ACG      | ACA      | T      | T      | synonymous     | 1      |
