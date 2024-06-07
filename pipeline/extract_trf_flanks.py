from collections import defaultdict
from pybedtools import BedTool
import sys
import argparse
import os

def get_scaffold_sizes(asm_fa):
    asm_fai = asm_fa + '.fai'
    if not os.path.exists(asm_fai):
        sys.exit('{} does not exist'.format(asm_fai))

    sizes = {}
    with open(asm_fai, 'r') as ff:
        for line in ff:
            scaffold, size = line.rstrip().split()[:2]
            sizes[scaffold] = int(size)
    return sizes

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("trf_tsv", type=str, help="screened trf tsv")
    parser.add_argument("asm_fa", type=str, help="indexed assembly fasta")     
    parser.add_argument("out_fa", type=str, help="out_fa")
    parser.add_argument("--flank_size", type=int, default=500, help="flank size. Default:500")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    bed_str = ''
    scaffold_sizes = get_scaffold_sizes(args.asm_fa)

    asm_loci = set()
    with open(args.trf_tsv, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split('\t')
            seq = cols[0]
            hit = cols[1:]
            asm_locus = '_'.join((seq, str(hit[0]), str(hit[1])))
            if asm_locus in asm_loci:
                continue
            asm_loci.add(asm_locus)
            start = int(hit[0]) - args.flank_size
            end = int(hit[1]) + args.flank_size

            #if start < 0 or end > scaffold_sizes[seq]:
            #    continue
            left = str(start - 1), str(int(hit[0]) - 1)
            right = hit[1], str(end)
            if int(left[0]) < 0 or int(right[1]) > scaffold_sizes[seq]:
                #print('overboard', asm_locus, left, right)
                continue
            bed_str += '{}\n'.format('\t'.join((seq, left[0], left[1], asm_locus + '_L')))
            bed_str += '{}\n'.format('\t'.join((seq, right[0], right[1], asm_locus + '_R')))

    flanks = BedTool(bed_str, from_string=True).sort()
    flanks = flanks.sequence(fi=args.asm_fa, name=True)

    with open(args.out_fa, 'w') as out:
        out.write('{}'.format(open(flanks.seqfn).read()))

main()
