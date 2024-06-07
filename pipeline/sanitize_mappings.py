from pybedtools import BedTool
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("tr_tsv", type=str, help="tr tsv")
    parser.add_argument("asm_mappings", type=str, help="assembly mappings tsv")
    parser.add_argument("out", type=str, help="output tsv")
    parser.add_argument("--log", type=str, help="log")
    parser.add_argument("--label", type=str, help="label, for appending haplotype to each row")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    tr_bed = BedTool(args.tr_tsv)
    asm_bed = BedTool(args.asm_mappings)
 
    tr_ncols = 8
    asm_ncols = 10
    loci_matched = set()
    loci_unmatched = set()
    loci_redundant = set()
    loci_initial = set()
    with open(args.out, 'w') as out:
        for cols in tr_bed.intersect(asm_bed, f=1.0, wo=True):
            tr_cols = cols[:8]
            asm_cols = cols[8:-1]
            tr_scaffold = tr_cols[4].split(':')[0]
            asm_scaffold = asm_cols[3]
            locus = tuple(tr_cols[:3])
            loci_initial.add(locus)

            if tr_scaffold == asm_scaffold:
                if not locus in loci_matched:
                    if args.label:
                        tr_cols.append(args.label)
                    out.write('{}\n'.format('\t'.join(tr_cols)))
                    loci_matched.add(locus)
                else:
                    loci_redundant.add(locus)
            else:
                loci_unmatched.add(locus)

    if args.log:
        with open(args.log, 'w') as out:
            out.write('loci initial:{}\n'.format(len(loci_initial)))
            out.write('loci matched:{}\n'.format(len(loci_matched)))
            out.write('loci unmatched:{}\n'.format(len(loci_unmatched)))
            out.write('loci multimapped:{}\n'.format(len(loci_unmatched & loci_matched)))
            out.write('loci redundant:{}\n'.format(len(loci_redundant)))

main()
