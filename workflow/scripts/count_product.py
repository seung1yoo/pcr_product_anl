

from Bio import SeqIO

def main(args):

    seq_dic = dict()
    for record in SeqIO.parse(open(args.ifn), 'fasta'):
        seq_dic.setdefault(str(record.seq), 0)
        seq_dic[str(record.seq)] += 1

    outfh = open(args.ofn, 'w')
    header = ['product_seq', 'count']
    outfh.write('{0}\n'.format('\t'.join(header)))
    for seq, cnt in seq_dic.items():
        outfh.write(f'{seq}\t{cnt}\n')
    outfh.close()



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifn')
    parser.add_argument('--ofn')
    args = parser.parse_args()
    main(args)
