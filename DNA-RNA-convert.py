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

# function to ranslate DNA sequence to mRNA to protein
def dna_to_protein(dna_sequence):
    # replacing 'T' with 'U' to create mRNA sequence
    mrna_sequence = dna_sequence.replace('T', 'U')
    
    # translate mRNA to protein
    protein = []
    for i in range(0, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]
        amino_acid = codon_table.get(codon, '')  # Get corresponding amino acid
        if amino_acid == 'Stop':
            break
        protein.append(amino_acid)
    
    return ''.join(protein)

# usage:
dna_sequence = input("Enter a DNA sequence (multiple of 3 nucleotides and contain only letters of A, T, C, G) Ex: ATGTTCGGA\n\n\n ").upper()

# validate (length multiple 3 and only letters A, T, C, G)
if len(dna_sequence) % 3 == 0 and all(base in 'ATCG' for base in dna_sequence):
    protein = dna_to_protein(dna_sequence)
    print(f"\nProtein sequence: {protein}\n")
else:
    print("\nInvalid DNA sequence, must be multiple of 3 and contains only A, T, C, G.\n")


