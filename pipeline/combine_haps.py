import argparse
import os
from glob import glob
from pybedtools import BedTool
import pysam
from collections import defaultdict, Counter
import sys
import pandas as pd
from graph import Graph
from operator import itemgetter
import re
import multiprocessing as mp
import gzip

def pick_coord(loci):
    coords = []
    for sample in loci.keys():
        for locus in loci[sample]:
            coords.append(tuple(locus[:3]))
    counts = Counter(coords)
    return counts.most_common(1)[0][0]

def pick_motif(loci):
    motifs = []
    #all_motifs = []
    for sample in loci.keys():
        for locus in loci[sample]:
            motifs.append(locus[3])

    counts = Counter(motifs)
    freqs = ';'.join(['{}({})'.format(count[0], count[1]) for count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))])
    return counts.most_common(1)[0][0], motifs

def is_same_repeat(rep1, rep2):
    if len(rep1) != len(rep2):
        return False
    if rep1 == rep2:
        return True

    for i in range(0, len(rep1), 1):
        rep = rep1[i:] + rep1[:i]
        if rep == rep2:
            return True

    return False

def count_motifs(motifs):
    counts = {}
    alias = {}
    used = set()
    for i in range(len(motifs)-1):
        if i in used:
            continue
        counts[i] = 1
        for j in range(i + 1, len(motifs)):
            if j in used:
                continue
            if is_same_repeat(motifs[i], motifs[j]):
                alias[motifs[j]] = motifs[i]
                used.add(j)
                counts[i] += 1
    return Counter({motifs[i]:counts[i] for i in counts}), alias

def update_motifs(loci):
    for locus in loci:
        all_motifs = sorted([m for m in locus[6].split(';') if m != '-'])

        counts, alias = count_motifs(all_motifs)
        consensus_motif = counts.most_common(1)[0][0]

        freqs = ';'.join(['{}({})'.format(count[0], count[1]) for count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))])
        locus[3] = consensus_motif
        locus.append(freqs)

        # update individual motifs
        motifs_new = []
        for motif in locus[6].split(';'):
            if motif in alias:
                motifs_new.append(alias[motif])
            else:
                motifs_new.append(motif)
        locus[6] = ';'.join(motifs_new)

def pick_locus(loci, multi_locus_samples, consensus_coord):
    consensus_coords = [int(consensus_coord[1]), int(consensus_coord[2])]
    consensus_span = consensus_coords[1] - consensus_coords[0] + 1
 
    for sample in multi_locus_samples:
        locus_olaps = []
        for locus in loci[sample]:
            olap = 0
            coord = int(locus[1]), int(locus[2])
            if coord == consensus_coords:
                # bonus 1.0 score to perfect match
                olap = 2.0
            else:
                olap = (coord[1] - coord[0] + 1) / consensus_span
            locus_olaps.append((locus, olap))

        if not locus_olaps:
            continue
        locus_olaps_sorted = sorted(locus_olaps, key=itemgetter(1), reverse=True)
        loci[sample] = [locus_olaps_sorted[0][0]]
        
def get_loci_across_samples(ccs, index_to_locus, samples):
    all_loci = []

    for cc in ccs:
        loci = defaultdict(list)
        for index in cc:
            locus = index_to_locus[index][:-1]
            sample = index_to_locus[index][-1]
            loci[sample].append(locus)
        consensus_coord = pick_coord(loci)
        consensus_motif, all_motifs = pick_motif(loci)

        # pick locus if sample has more than 1 locus
        multi_locus_samples = [s for s in loci.keys() if len(loci[s]) > 1]
        if multi_locus_samples:
            pick_locus(loci, multi_locus_samples, consensus_coord)

        locus = list(consensus_coord)
        locus.append(consensus_motif)
        sizes = []
        copy_nums = []
        motifs = []
        for sample in samples:
            if sample in loci and loci[sample]:
                sizes.append(loci[sample][0][5])
                copy_nums.append(loci[sample][0][6])
                motifs.append(loci[sample][0][3])
            else:
                sizes.append('-')
                copy_nums.append('-')
                motifs.append('-')
        locus.append(';'.join(copy_nums))
        locus.append(';'.join(sizes))
        locus.append(';'.join(motifs))

        all_loci.append(locus)
    
    return all_loci
        
def find_ccs(olaps):
    all_loci = set([o[0] for o in olaps] + [o[1] for o in olaps])
    print('nloci', len([o[0] for o in olaps] + [o[1] for o in olaps]), len(all_loci))
    
    g = Graph(len(all_loci))
    i = 0
    used = {}
    for olap in olaps:
        index = []
        for locus in olap:
            if locus in used:
                index.append(used[locus])
            else:
                index.append(i)
                used[locus] = i
                i += 1
        g.addEdge(index[0], index[1])
    
    ccs = g.connectedComponents()

    index_to_locus = {}
    for locus, index in used.items():
        index_to_locus[index] = locus

    return ccs, index_to_locus

def split_tasks(args, n):
    k, m = divmod(len(args), n)
    batches = (args[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
    return [b for b in batches if b]

def overlap(samples, regions, nprocs, outdir=None):
    results = []
    if nprocs > 1:
        tasks = []
        for region in regions:
            tasks.append((samples, region, outdir))
        batches = list(split_tasks(tasks, nprocs))
        pool = mp.Pool(processes=nprocs)
        batch_results = pool.map(merge_bed_worker, batches)
        results = [result for batch in batch_results for result in batch]
    else:
        for region in regions:
            loci = merge_samples(samples, region)
            results.append((loci, region))

    loci = list(sum([tuple(r[0]) for r in results], ()))

    return loci

def merge_bed_worker(tasks):
    results = []
    for samples, region, outdir in tasks:
        sample_names = [s[0] for s in samples]
        loci = merge_samples(samples, region)
        post_process(loci, sample_names)
        results.append((loci, region))
        if outdir is not None:
            out_tsv = '{}/{}.tsv'.format(outdir, region)
            output(loci, out_tsv, sample_names)

    return results

def merge_samples(samples, region, ncols=9, min_f=0.7):
    beds = []
    samples_sorted = sorted(samples, key=itemgetter(0))
    used = defaultdict(set)
    for sample, tsv in samples_sorted:
        if region:
            bed = BedTool(tsv).tabix_intervals(region)
        else:
            bed = BedTool(tsv)
        beds.append(bed)

    all_olaps = []
    for i in range(len(beds)-1):
        sample = samples[i][0]
        if not sample in used:
            a = beds[i]
        else:
            bed_str = ''
            for locus in used[sample]:
                bed_str += '{}\n'.format('\t'.join(locus))
            a = beds[i].subtract(BedTool(bed_str, from_string=True).sort(), f=1.0, r=True)

        b = [beds[j].fn for j in range(i+1, len(beds))]
        olaps = merge_bed(a, b)

        for olap in olaps:
            for locus in olap:
                used[locus[-1]].add(locus)
        
        all_olaps.extend(olaps)

    ccs, index_to_locus = find_ccs(all_olaps)

    loci = get_loci_across_samples(ccs, index_to_locus, [s[0] for s in samples_sorted])

    return loci

def merge_bed(bed1, bed2, max_n=510, ncols=9, min_f=0.7):
    olaps = []
    for b2 in [bed2[i : i + max_n] for i in range(0, len(bed2), max_n)]:
        for cols in bed1.intersect(b2, f=min_f, sorted=True, r=True, wo=True):
            locus1 = tuple(cols[:ncols])
            # 1 extra column in between 2 loci when b is multiple bed files
            if len(b2) > 1:
                locus2 = tuple(cols[ncols+1:-1])
            else:
                locus2 = tuple(cols[ncols:-1])
            olaps.append((locus1, locus2))
    
    return olaps

def add_size_changes(loci):
    for locus in loci:
        ref_size = int(locus[2]) - int(locus[1]) + 1

        all_changes = []
        for sizes in locus[5].split(';'):
            changes = []
            for size in sizes.split(','):
                if size != '-':
                    changes.append(int(size) - ref_size)
                else:
                    changes.append('-')
            all_changes.append(','.join(list(map(str, changes))))
        locus.append(';'.join(all_changes))

def add_max_size_change(loci):
    for locus in loci:
        ref_size = int(locus[2]) - int(locus[1]) + 1

        changes = []
        for sizes in locus[5].split(';'):
            for size in sizes.split(','):
                if size != '-':
                    changes.append(abs(float(size) - ref_size))
        if not changes:
            max_change = 0
        else:
            max_change = max(changes)
        # keep using int for df deletion
        locus.append(max_change)

def add_genotype_counts(loci, samples):
    for locus in loci:
        num_haps_genotyped = 0
        uniq_samples = set()
        for sizes, sample in zip(locus[5].split(';'), samples):
            for size in sizes.split(','):
                if not '-' in size:
                    num_haps_genotyped += 1
                    uniq_samples.add(sample.split('.')[0])
        locus.append(str(len(uniq_samples)))
        locus.append(str(num_haps_genotyped))

def update_copy_nums(loci):
    for locus in loci:
        motif_size = len(locus[3])

        all_copy_nums = []
        for sizes in locus[5].split(';'):
            copy_nums = []
            for size in sizes.split(','):
                if size != '-':
                    copy_nums.append('{:.1f}'.format(float(size)/motif_size))
                else:
                    copy_nums.append('-')
            all_copy_nums.append(','.join(list(map(str, copy_nums))))
        locus[4] = ';'.join(all_copy_nums)

def post_process(loci, samples):
    print('updating copy numbers...')
    update_copy_nums(loci)
    print('adding size changes...')
    add_max_size_change(loci)
    print('adding genotype_counts...')
    add_genotype_counts(loci, samples)
    print('updating motifs...')
    update_motifs(loci)

def output(loci, out_file, samples, remove_unchanged=False, hide_motifs=False):
    out = open(out_file, 'w')

    calls = ';'.join(samples)
    nsamples = len(set([s.split('.')[0] for s in samples]))
    out.write('# {} {} {}\n'.format(calls, nsamples, len(samples)))
    header = ('chrom', 'start', 'end', 'motif', 'copy_numbers', 'sizes', 'motifs', 'max_change', 'num_samples', 'num_calls', 'motif_frequency')

    data = []
    for locus in loci:
        cols = [locus[0], int(locus[1]), int(locus[2])] + locus[3:]
        data.append(cols)
    df = pd.DataFrame(data, columns=header)
    df.drop_duplicates(inplace=True)
    if hide_motifs:
        df.drop(columns='motifs', inplace=True)
    if remove_unchanged:
        df = df[df.max_change != 0]
    df.sort_values(by=['chrom', 'start', 'end'], inplace=True)
    df.to_csv(out, sep='\t', index=False)

def get_regions(tsv_gz, skip_chroms):
    bands = []
    with gzip.open(tsv_gz, 'rt') as ff:
        for line in ff:
            if line[0] == '#':
                continue
            cols = line.rstrip().split()
            if '_' in cols[0] or cols[0] in skip_chroms or cols[-1] == 'acen':
                continue
            bands.append('{}:{}-{}'.format(cols[0], cols[1], cols[2]))

    return bands

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("topdir", type=str, help="top directory where all data reside")
    parser.add_argument("suffix", type=str, help="suffix of parental results following sample name e.g HG002.suffix")
    parser.add_argument("out", type=str, help="output file")
    parser.add_argument("--ideo_gz", type=str, help="gzipped ideogram file")
    parser.add_argument("--nprocs", type=int, default=1, help="number of processes (Default=1)")
    parser.add_argument("--skip_chroms", type=str, default="chrY,chrM", help="comma-separated list of skip chromosomes, e.g. chrY,chrM")
    parser.add_argument("--samples", type=str, help="comma-separated list of samples")
    parser.add_argument("--regions", type=str, help="comma-separated list of UCSC-format regions")
    parser.add_argument("--hide_motifs", action='store_true', help="hide motif details")
    parser.add_argument("--outdir", type=str, help="output directory for intermediate tsvs") 
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    skip_chroms = args.skip_chroms.split(',')
    samples_only = []
    if args.samples:
        samples_only = args.samples.split(',')

    all_bed_lines = []
    all_samples = []
    sample_dirs = {}
    for d in glob('{}/*'.format(args.topdir)):
        sample = os.path.basename(d)
        if samples_only and not sample in samples_only:
            continue
        sample_dirs[sample] = d
        tsvs = glob('{}/{}.*.{}'.format(d, sample, args.suffix))
        for tsv in tsvs:
            prefix = '.'.join(os.path.basename(tsv).split('.')[:2])
            all_samples.append((prefix, tsv))
    
    if all_samples:
        all_samples.sort(key=itemgetter(0))
        regions = []
        if args.regions:
            regions = args.regions.split(',')
        elif args.ideo_gz:
            regions = get_regions(args.ideo_gz, skip_chroms)

        loci = overlap(all_samples, regions, args.nprocs, outdir=args.outdir)
        samples = [s[0] for s in all_samples]

        output(loci, args.out, samples, hide_motifs=args.hide_motifs)

main()
