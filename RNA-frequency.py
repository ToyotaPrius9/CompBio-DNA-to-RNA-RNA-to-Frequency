# dictionary for RNA codons & amino acids
codon_table = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# reversing codon table to map amino acids-codons
amino_acid_to_codons = {}
for codon, amino_acid in codon_table.items():
    if amino_acid not in amino_acid_to_codons:
        amino_acid_to_codons[amino_acid] = []
    amino_acid_to_codons[amino_acid].append(codon)

# function for finding RNA codon frequencies for specific amino acids in DNA sequence
def codon_frequency(dna_sequence, amino_acids):
    # replacing 'T' with 'U' for making mRNA sequence
    mrna_sequence = dna_sequence.replace('T', 'U')
    
    # dictionary to store codon frequencies
    codon_count = {amino_acid: {codon: 0 for codon in amino_acid_to_codons[amino_acid]} for amino_acid in amino_acids}
    
    # loop through mRNA sequence in codons (triplet)
    for i in range(0, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]
        # finding the amino acid for codon
        for amino_acid in amino_acids:
            if codon in amino_acid_to_codons[amino_acid]:
                codon_count[amino_acid][codon] += 1
                
    return codon_count




# usage
dna_sequence = input("\nEnter a DNA sequence (example: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG):\n\n ").upper()
amino_acids = input("\nenter up to 3 amino acids (single-letter codes, like A, F, L):\n\n ").upper().split()

# validate
valid_amino_acids = set(codon_table.values())

# usage
amino_acids = [aa for aa in amino_acids if aa in valid_amino_acids][:3]  # limit to 3 aa
if amino_acids:
    # get freq
    codon_count = codon_frequency(dna_sequence, amino_acids)

    # fancy print
    for amino_acid, codons in codon_count.items():
        print(f"\nAmino acid: {amino_acid}")
        for codon, count in codons.items():
            print(f"  Codon {codon}: {count}")
else:
    print("\namino acids invalid\n\n")
