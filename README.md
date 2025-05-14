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

## Codon Comparison
