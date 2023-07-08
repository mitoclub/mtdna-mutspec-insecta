from Bio import SeqIO

# Separate IDs to single files (extracting only protein sequence)

fin = "../../raw/blattodea_genbank/cockroach_complete_refseq.gb"
#fin = "../../raw/blattodea_genbank/TermiteRefSeqs.gb"


for entry in SeqIO.parse(fin, 'gb'):
    for feat in entry.features:
        if feat.type == 'CDS':
            if 'gene' in feat.qualifiers:
                with open(f"../../interim/ForDolphin/separate_ids_cock/{entry.id}_{feat.qualifiers['gene'][0]}.faa", "w") as fout:
                    fout.write(f">{entry.id} [{entry.description}]\n{feat.qualifiers['translation'][0]}")

