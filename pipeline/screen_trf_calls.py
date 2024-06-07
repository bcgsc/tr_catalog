import argparse
import sys
import re
from pybedtools import BedTool

def parse_trf(trf_output, max_motif_len):
    bed_str = ''
    with open(trf_output, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split()
            if not cols:
                continue
            if cols[0] == 'Sequence:':
                seq = cols[1]
                m = re.search(r"(.+):(\d+)", seq)
                offset = 0
                if m:
                    seq = m.group(1)
                    offset = int(m.group(2))
            elif len(cols) == 15:
                if len(set(cols[13])) > 1 and len(cols[13]) <= max_motif_len:
                    cols[0] = str(int(cols[0]) + offset)
                    cols[1] = str(int(cols[1]) + offset)
                    bed_str += '{}\n'.format('\t'.join([seq] + cols))

    return BedTool(bed_str, from_string=True).sort()

def screen_trf_results(bed, out_file, ncols=16, f=0.8):
    keep = set()
    remove = set()
    for cols in bed.intersect(bed, sorted=True, f=f, r=True, wo=True):
        cols1 = tuple(cols[:ncols])
        cols2 = tuple(cols[ncols:ncols*2])

        if cols1 == cols2:
            keep.add(cols1)
        elif cols1[0] == cols2[0] and int(cols1[1]) <= int(cols2[1]):
            #print('bb', cols1, cols2)
            if int(cols1[3]) < int(cols2[3]):
                remove.add(cols2)
            elif int(cols2[3]) < int(cols1[3]):
                remove.add(cols1)
            elif int(cols1[8]) > int(cols2[8]):
                remove.add(cols2)
            elif int(cols2[8]) > int(cols1[8]):
                remove.add(cols1)

    #print('aa', len(keep), len(remove), len(keep-remove))
    bed_str = ''
    for cols in keep - remove:
        bed_str += '{}\n'.format('\t'.join(cols))

    BedTool(bed_str, from_string=True).sort().moveto(out_file)

def output(results, out_file):
    with open(out_file, 'w') as out:
        for seq in sorted(results.keys()):
            for cols in results[seq]:
                out.write('{}\n'.format('\t'.join([seq] + list(map(str, cols)))))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("trf", type=str, help="input trf")
    parser.add_argument("out", type=str, help="output")
    parser.add_argument("--max_motif_len", type=int, help="maximum motif len", default=100)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    bed = parse_trf(args.trf, args.max_motif_len)
    screen_trf_results(bed, args.out)

main()
