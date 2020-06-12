import sys
import re
import pandas as pd
from itertools import compress

def find_KK_seq_from_fasta(fasta_file, output_file):
    dilysine = r'(KK...|K.K..)'
    proteins = {}
    with open(fasta_file) as proteins_fasta:
        for line in proteins_fasta:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                protein_name = line[1:]
                if protein_name not in proteins:
                    proteins[protein_name] = []
                continue
            sequence: str = line
            proteins[protein_name].append(sequence)

    for key, value in proteins.items():
        proteins[key] = "".join(value).upper()
    proteins_df = pd.DataFrame(proteins.items(), columns=['protein', 'aa_sequence'])
    proteins_df['seq_length'] = proteins_df.aa_sequence.apply(len)

    matches = []
    positions = []
    for seq in proteins.values():
        match = re.findall(dilysine, seq)
        pos = [seq.find(pat) + 1 for pat in match]
        mask = [len(seq) - p == 4 for p in pos]
        match = list(compress(match, mask))
        pos = list(compress(pos, mask))
        matches.append(match)
        positions.append(pos)
    proteins_df['KK_seq'] = matches
    proteins_df['KK_positions'] = positions

    proteins_df.to_csv(output_file)
    return (proteins_df)