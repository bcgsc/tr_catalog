import pysam
import sys
from collections import defaultdict
import argparse
from pybedtools import BedTool
import os

def type_trf_cols(cols):
    return list(map(int, cols[:3])) + [float(cols[3])] + list(map(int, cols[4:12])) + [float(cols[12])] + cols[13:]

def parse_trf(trf_output):
    results = defaultdict(list)
    with open(trf_output, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split()
            if not cols:
                continue
            if cols[0] == 'Sequence:':
                seq = cols[1]
            elif len(cols) == 15:
                if len(cols[13]) > 1:
                    results[seq].append(type_trf_cols(cols))

    return results

def parse_trf_tsv(trf_output):
    results = defaultdict(list)
    with open(trf_output, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split('\t')
            results[cols[0]].append(cols[1:])

    return results

def get_motifs(trf_out):
    motifs = {}
    if os.path.splitext(trf_out)[1] == '.dat':
        trf_results = parse_trf(trf_out)
    elif os.path.splitext(trf_out)[1] == '.tsv':
        trf_results = parse_trf_tsv(trf_out)
    for seq, hits in trf_results.items():
        for hit in hits:
            asm_locus = '_'.join((seq, str(hit[0]), str(hit[1])))
            if not asm_locus in motifs or len(hit[-2]) < len(motifs[asm_locus]):
                motifs[asm_locus] = hit[-2]

    return motifs

def good_end(end_tuple):
    return end_tuple[0] == 0 or (end_tuple[0] >= 4 and end_tuple[0] <= 5 and end_tuple[1] <= 10)

def get_span(aln):
    flank = aln.query_name.split(':')[0].rsplit('_', 1)[1]                                                
    if flank == 'L':
        end_tuple = aln.cigartuples[0] if not aln.is_reverse else aln.cigartuples[-1]
    else:
        end_tuple = aln.cigartuples[-1] if not aln.is_reverse else aln.cigartuples[0]

    if not good_end(end_tuple):
        return None
    
    if flank == 'L':
        if not aln.is_reverse:
            start, end = aln.reference_end + 1, None
        else:
            start, end = None, aln.reference_start
    elif flank == 'R':
        if not aln.is_reverse:
            start, end = None, aln.reference_start
        else:
            start, end = aln.reference_end + 1, None
    else:
        return None

    return (aln.reference_name, start, end)

def are_different(alns):
    by_flank = defaultdict(list)
    for aln in alns:
        flank = aln.query_name.rsplit('_', 1)[1]
        by_flank[flank].append(aln)

    if len(by_flank['L']) == len(by_flank['R']):
        uniq = True
        for flank in ('L', 'R'):
            if len(set([aln.to_string() for aln in by_flank[flank]])) != 1:
                uniq = False

        if uniq:
            return (by_flank['L'][0], by_flank['R'][0])

    return False

def get_mappings(bam):
    loci = defaultdict(list)
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_supplementary:
            continue

        locus = aln.query_name.split('::')[0].rsplit('_', 1)[0]
        loci[locus].append(aln)
    
    mappings = defaultdict(list)
    failed = {}
    for locus, alns in loci.items():
        if len(alns) == 1:
            print('failed_only_1_aligned', locus, len(alns))
            failed[locus] = 'only_1_aligned'
            continue

        aln1, aln2 = alns[:2]
        if len(alns) > 2:
            diff_alns = are_different(alns)
            if diff_alns:
                aln1, aln2 = diff_alns
            else:
                print('failed_same_aligns', locus, len(alns))
                failed[locus] = 'same_aligns'
                continue

        if aln1.reference_name == aln2.reference_name and aln1.is_reverse != aln2.is_reverse:
            print('failed_same_strand', locus, aln1.reference_name, aln1.is_reverse, aln2.reference_name, aln2.is_reverse)
            failed[locus] = 'same_strand'
            continue

        if aln1.query_name == aln2.query_name:
            print('failed_same_query', locus, aln1.query_name, aln1.reference_start, aln2.query_name, aln2.reference_start)
            failed[locus] = 'same_query'
            continue

        if aln1.reference_name != aln2.reference_name:
            print('diff chrom', locus, aln1.reference_name, aln2.reference_name)
            failed[locus] = 'diff_chrom'
            continue

        if aln1.has_tag('XA') or aln2.has_tag('XA'):
            print('warning_has_xa', locus, aln1.query_name, aln1.has_tag('XA'), aln2.query_name, aln2.has_tag('XA'))
            #continue

        asm_pos = list(map(int, locus.rsplit('_', 2)[-2:]))
        asm_span = asm_pos[1] - asm_pos[0] + 1

        pos1 = get_span(aln1)
        pos2 = get_span(aln2)
        if pos1 is None or pos2 is None:
            failed[locus] = 'no_span'
            if pos1 is None:
                print('failed_no_span', locus, aln1.reference_start, aln1.reference_end, aln1.cigarstring)
            else:
                print('failed_no_span', locus, aln2.reference_start, aln2.reference_end, aln2.cigarstring)
            continue
        strand = '-' if aln1.is_reverse else '+'

        chrom = pos1[0]
        start = None
        end = None
        if pos1[1] is None and pos2[2] is None:
            start, end = pos2[1], pos1[2]
        elif pos1[2] is None and pos2[1] is None:
            start, end = pos1[1], pos2[2]

        if start is not None and end is not None:
            if end > start:
                chrom_span = end - start + 1
                if chrom_span / asm_span > 10:
                    continue
                mappings[(chrom, start, end, chrom_span)].append((locus.rsplit('_', 2)[0], asm_pos[0], asm_pos[1], asm_span, strand))

            else:
                failed[locus] = 'start_after_end'
                print('failed start_after_end', locus, chrom, start, end, start-end)

    return mappings, failed

def reverse_complement(seq):
    """Reverse complements sequence string"""
    complement = str.maketrans("agtcAGTC", "tcagTCAG")
    return seq[::-1].translate(complement)

def output(mappings, motifs, out_file):
    bed_str = ''
    for ref_locus, asm_mappings in mappings.items():
        for asm_mapping in asm_mappings:
            motif = motifs['{}_{}_{}'.format(asm_mapping[0], asm_mapping[1], asm_mapping[2])]
            if asm_mapping[-1] == '-':
                motif = reverse_complement(motif)
            asm_locus = '{}:{}-{}:{}'.format(asm_mapping[0], asm_mapping[1], asm_mapping[2], asm_mapping[4])
            asm_span = asm_mapping[-2]
        
            cols = list(ref_locus[:3]) + [motif, asm_locus, asm_span, '{:.1f}'.format(asm_span / len(motif)), ref_locus[-1]]
            bed_str += '{}\n'.format('\t'.join(list(map(str, cols))))

    BedTool(bed_str, from_string=True).sort().saveas(out_file)

def output_fails(failed, out_bed):
    bed_str = ''
    for locus, reason in failed.items():
        cols = locus.split('_')
        cols.append(reason)
        bed_str += '{}\n'.format('\t'.join(cols))

    BedTool(bed_str, from_string=True).saveas(out_bed)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam")
    parser.add_argument("trf_out", type=str, help="trf output")
    parser.add_argument("out", type=str, help="out")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam)
    mappings, failed = get_mappings(bam)
    
    motifs = get_motifs(args.trf_out)

    output(mappings, motifs, args.out)
    prefix = os.path.splitext(args.out)[0]
    output_fails(failed, '{}.failed.bed'.format(prefix))

main()
