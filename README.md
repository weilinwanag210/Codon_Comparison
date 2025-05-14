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
One aa in reference will miss or be replaced by another amino acid or be inserted _XXXXX_in variant sequence with a fixed probablity
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

## Codon Comparison
