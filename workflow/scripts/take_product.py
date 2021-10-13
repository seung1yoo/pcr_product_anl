
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip

class Product:
    def __init__(self, args):
        self.primer5 = Seq(args.primer5)
        self.primer3 = Seq(args.primer3)
        self.infqgz = args.infqgz
        self.prefix = args.prefix
        #
        self.find_product(self.primer5, self.primer3, f"{self.prefix}.product")
        #
    def find_product(self, primer5, primer3, outprefix):
        outfna = open(f"{outprefix}.fna", 'w')
        outfaa = open(f"{outprefix}.faa", 'w')
        outtsv = open(f"{outprefix}.tsv", 'w')
        #
        for record in SeqIO.parse(gzip.open(self.infqgz, 'rt'), 'fastq'):
            loc5 = record.seq.find(str(primer5))
            loc3 = record.seq.find(str(primer3))
            if not loc5 in [-1] and not loc3 in [-1]:
                product_seq = str(record.seq[loc5+len(primer5):loc3])
                if not len(product_seq) in [33]:
                    continue
                if '*' in Seq(product_seq).translate():
                    continue
                #fna
                product_record = SeqRecord(
                    Seq(product_seq),
                    id=f"{record.id}_product",
                    description=f"strand=sense"
                )
                SeqIO.write(product_record, outfna, 'fasta')
                #faa
                product_record = SeqRecord(
                    Seq(product_seq).translate(),
                    id=f"{record.id}_product",
                    description=f"strand=sense"
                )
                SeqIO.write(product_record, outfaa, 'fasta')
                #
                items = [record.id]
                items.append('sense')
                items.append('found_both')
                items.append(loc5)
                items.append(loc3)
                outtsv.write('{0}\n'.format('\t'.join([str(x) for x in items])))
            elif not loc5 in [-1] and loc3 in [-1]:
                items = [record.id]
                items.append('sense')
                items.append('found_5primer')
                items.append(loc5)
                items.append(loc3)
                items.append(record.seq)
                outtsv.write('{0}\n'.format('\t'.join([str(x) for x in items])))
            elif loc5 in [-1] and not loc3 in [-1]:
                items = [record.id]
                items.append('sense')
                items.append('found_3primer')
                items.append(loc5)
                items.append(loc3)
                items.append(record.seq)
                outtsv.write('{0}\n'.format('\t'.join([str(x) for x in items])))
        #
        for record in SeqIO.parse(gzip.open(self.infqgz, 'rt'), 'fastq'):
            loc5 = record.seq.find(str(primer5.reverse_complement()))
            loc3 = record.seq.find(str(primer3.reverse_complement()))
            if not loc5 in [-1] and not loc3 in [-1]:
                product_seq = str(record.seq[loc3+len(primer3):loc5].reverse_complement())
                if not len(product_seq) in [33]:
                    continue
                if '*' in Seq(product_seq).translate():
                    continue
                #fna
                product_record = SeqRecord(
                    Seq(product_seq),
                    id=f"{record.id}_product",
                    description=f"strand=antisense"
                )
                SeqIO.write(product_record, outfna, 'fasta')
                #faa
                product_record = SeqRecord(
                    Seq(product_seq).translate(),
                    id=f"{record.id}_product",
                    description=f"strand=antisense"
                )
                SeqIO.write(product_record, outfaa, 'fasta')
                #
                items = [record.id]
                items.append('antisense')
                items.append('found_both')
                items.append(loc5)
                items.append(loc3)
                outtsv.write('{0}\n'.format('\t'.join([str(x) for x in items])))
            elif not loc5 in [-1] and loc3 in [-1]:
                items = [record.id]
                items.append('antisense')
                items.append('found_5primer')
                items.append(loc5)
                items.append(loc3)
                items.append(record.seq)
                outtsv.write('{0}\n'.format('\t'.join([str(x) for x in items])))
            elif loc5 in [-1] and not loc3 in [-1]:
                items = [record.id]
                items.append('antisense')
                items.append('found_3primer')
                items.append(loc5)
                items.append(loc3)
                items.append(record.seq)
                outtsv.write('{0}\n'.format('\t'.join([str(x) for x in items])))
        outfna.close()
        outfaa.close()
        outtsv.close()

def main(args):


    product = Product(args)



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--primer5')
    parser.add_argument('--primer3')
    parser.add_argument('--infqgz')
    parser.add_argument('--prefix')
    args = parser.parse_args()
    main(args)


