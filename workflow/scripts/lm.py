

import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm

from Bio import SeqIO

def run(fn, prefix, seqtype):
    df = pd.read_csv(fn, comment=">", names=['site'])
    count_df = lm.alignment_to_matrix(sequences=df['site'].values,
                                      to_type='counts', characters_to_ignore='*')
    count_df.to_csv(f"{prefix}.count.tsv", index=True, header=True, sep="\t")
    #print(count_df)
    prob_df = lm.transform_matrix(count_df, from_type='counts', to_type='probability')
    prob_df.to_csv(f"{prefix}.probability.tsv", index=True, header=True, sep="\t")
    #print(prob_df)
    if seqtype in ['fna']:
        logo = lm.Logo(prob_df, color_scheme="colorblind_safe", figsize=(10,3))
        #logo.ax.set_xlabel("")
        #logo.ax.set_ylabel("")
    elif seqtype in ['faa']:
        logo = lm.Logo(prob_df, color_scheme="dmslogo_funcgroup", figsize=(10,3))
    plt.savefig(f"{prefix}.png")

    outfh = open(f"{prefix}", "w")
    for record in SeqIO.parse(open(fn), 'fasta'):
        SeqIO.write(record, outfh, 'fasta')
    outfh.close()

def main(args):
    run(args.infa, f"{args.prefix}", args.seqtype)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--infa')
    parser.add_argument('--prefix')
    parser.add_argument('--seqtype', choices=['fna','faa'])
    args = parser.parse_args()
    main(args)
