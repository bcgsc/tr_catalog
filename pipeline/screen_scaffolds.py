import argparse
from operator import itemgetter
from collections import defaultdict
from pybedtools import BedTool

def parse_paf(paf):
    by_scaffolds = defaultdict(list)
    with open(paf, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split('\t')
            #qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, res_matches, alen, mapq = cols[:12]
            qname = cols[0]
            strand = cols[4]
            tname = cols[5]
            qlen, qstart, qend, tlen, tstart, tend, res_matches, alen, mapq = list(map(int, cols[1:4] + cols[6:12]))
            if mapq >= 60 and 'MT' not in qname and tname != "chrY":
                by_scaffolds[qname].append((qlen, qstart, qend, strand, tname, tlen, tstart, tend, res_matches, alen, mapq))

    return by_scaffolds

def screen_mappings(mappings, cen_bed=None):
    bed_str = ''
    ncols = None
    for scaffold in mappings:
        mappings_sorted = sorted(mappings[scaffold], key=itemgetter(9), reverse=True)
        mappings_kept = []
        chrom = None
        for i in range(len(mappings_sorted)):
            m = mappings_sorted[i]
            if i == 0:
                mappings_kept.append(m)
                chrom = m[4]
            else:
                if m[4] == chrom:
                    mappings_kept.append(m)
        
        for m in mappings_kept:
            cols = [m[4], m[6], m[7], scaffold, m[0], m[1], m[2], m[3], m[9], m[10]]
            bed_str += '{}\n'.format('\t'.join(list(map(str, cols))))
            if ncols is None:
                ncols = len(cols)

    bed = BedTool(bed_str, from_string=True).sort()
    if cen_bed is not None:
        bed.intersect(cen_bed, f=0.8, wo=True)
        bed = bed.intersect(cen_bed, f=0.8, v=True)
    #bed.saveas('bb.bed')

    bed_mm = find_multimaps(bed, ncols)
    #bed_mm.saveas('cc.bed')

    if bed_mm:
        bed = bed.subtract(bed_mm)
    #bed.saveas('dd.bed')

    return bed

def output_trf_inputs(bed, asm_fa, outdir):
    for cols in bed:
        bed = BedTool('{}\t{}\t{}'.format(cols[3], cols[5], cols[6]), from_string=True)
        bed = bed.sequence(fi=asm_fa)
        out_fa = '{}/{}:{}-{}.fa'.format(outdir, cols[3], cols[5], cols[6])
        with open(out_fa, 'w') as out:
            out.write('{}'.format(open(bed.seqfn).read()))

def find_multimaps(bed, ncols):
    mms = set()
    for cols in bed.intersect(bed, f=0.8, r=True, wo=True):
        segment1 = tuple(cols[:ncols])
        segment2 = tuple(cols[ncols:ncols*2])
        if segment1 == segment2:
            continue
        mms.add(segment1)
        mms.add(segment2)

    bed_str = ''
    for cols in mms:
        bed_str += '{}\n'.format('\t'.join(cols))
    return BedTool(bed_str, from_string=True)
 
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("paf", type=str, help="paf") 
    parser.add_argument("out", type=str, help="output mappings tsv")
    parser.add_argument("--asm", type=str, help="indexed assembly fasta")
    parser.add_argument("--trf_inputs", type=str, help="trf inputs dir")
    parser.add_argument("--cen", type=str)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    bed = screen_mappings(parse_paf(args.paf), cen_bed=args.cen)
    if args.asm and args.trf_inputs:
        output_trf_inputs(bed, args.asm, args.trf_inputs)
    bed.moveto(args.out)

main()
