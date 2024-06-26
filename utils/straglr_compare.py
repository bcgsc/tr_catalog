#!/usr/bin/env python
import argparse
from pybedtools import BedTool, create_interval_from_list
from collections import defaultdict, Counter
from scipy import stats
import numpy as np
from operator import itemgetter
import sys
import os
import re
import warnings
# turn off RuntimeWarning from stats.ttest_ind()
warnings.filterwarnings('ignore')
import pandas as pd

def create_straglr_bed(alleles):
    bed_str = ''
    loci = defaultdict(dict)
    max_num_alleles = 0
    for locus in alleles:
        if len(alleles[locus]) > max_num_alleles:
            max_num_alleles = len(alleles[locus])
        for allele in alleles[locus]:
            loci[locus][float(allele)] = ','.join(alleles[locus][allele])

    ncols = 5 + 2 * max_num_alleles
    for locus in loci:
        cols = list(locus)
        for allele in sorted(loci[locus]):
            cols.extend([str(allele), loci[locus][allele]])
        if len(cols) < ncols:
            cols.extend(['-'] * (ncols - len(cols)))
        bed_str += '{}\n'.format('\t'.join(cols))

    return BedTool(bed_str, from_string=True).sort()

def create_simple_bed(loci):
    bed_str = ''
    for locus in loci:
        bed_str += '{}\n'.format('\t'.join(list(map(str, locus[:3]))))

    return BedTool(bed_str, from_string=True).sort()

def parse_straglr_tsv(tsv, use_size=True, skip_chroms=None, old_version=False):
    alleles = {}
    with open(tsv, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                continue

            cols = line.rstrip().split('\t')
            if not old_version and cols[-1] != 'full':
                continue
            
            if skip_chroms is not None and cols[0] in skip_chroms:
                continue
            locus = cols[:4]
            # coordinate error in straglr output
            if int(locus[1]) > int(locus[2]):
                continue
            ref_size = int(cols[2]) - int(cols[1]) + 1
            if not use_size:
                ref_size = '{:.1f}'.format(ref_size / len(cols[3]))
            locus.append(str(ref_size))
            locus = tuple(locus)
            if not locus in alleles:
                alleles[locus] = defaultdict(list)
            allele = cols[-1] if old_version else cols[-2]
            if allele == '-' or allele == 'NA':
                continue
            if not use_size:
                allele = '{:.1f}'.format(float(allele) / len(cols[3]))
            
            size = cols[-5]
            copy_num = cols[-6]
            if old_version:
                size = cols[-3]
                copy_num = cols[-4]
            
            if use_size:
                alleles[locus][allele].append(size)
            else:
                alleles[locus][allele].append(copy_num)

    return create_straglr_bed(alleles)

def pick_best_olap(olaps):
    ''' pick best control that has the largest overlap with the test locus '''
    for loc in olaps:
        if len(olaps[loc]) > 1:
            fs = []
            for i in range(len(olaps[loc])):
                cols = olaps[loc][i]
                span1 = int(cols[2]) - int(cols[1]) + 1
                span2 = int(cols[11]) - int(cols[10]) + 1
                olap = min(float(cols[2]), float(cols[11])) - max(float(cols[1]), float(cols[10])) + 1
                f1 = olap/span1
                f2 = olap/span2
                f = (f1 + f2) / 2
                fs.append((i, f))
            fs_sorted = sorted(fs, key=itemgetter(1), reverse=True)
            best = fs_sorted[0][0]
            for i in range(len(olaps[loc])-1,-1,-1):
                if i != best:
                    del(olaps[loc][i])

def get_larger_alleles(ns, min_count=2):
    ''' pick the larger alleles from the heterzygous allele '''
    ns_sorted = sorted(ns)
    max_jump = 0
    index = None
    for i in range(len(ns)-1):
        jump = ns_sorted[i+1] - ns_sorted[i]
        if jump > max_jump:
            max_jump = jump
            index = i + 1

    larger_alleles = ns_sorted[index:]

    if len(larger_alleles) >= min_count:
        return larger_alleles
    else:
        return ns_sorted

def vs_each_control(test_bed, control_bed, pval_cutoff, min_expansion=0, min_support=0, label=None, from_catalog=False, use_size=False):
    expanded_loci = {}
    new_loci = []
    common_loci = []
    for cols in test_bed.intersect(control_bed, loj=True, wao=True, f=0.9).saveas('a1.bed'):
        if cols[9] == '.':
            new_loci.append(cols)
        else:
            common_loci.append(tuple(cols))
    for cols in test_bed.intersect(control_bed, loj=True, wao=True, F=0.9).saveas('a2.bed'):
        if cols[9] != '.':
            if not tuple(cols) in common_loci:
                common_loci.append(tuple(cols))
            if cols in new_loci:
                new_loci.remove(cols)

    for cols in new_loci:
        locus = (cols[0], int(cols[1]), int(cols[2]), cols[3], float(cols[4]))
        ref_size = locus[-1]
        test_alleles = []
        supports = {}

        for i in range(5, 8, 2):
            if cols[i] == '-' or float(cols[i]) - ref_size < min_expansion:
                continue
            test_calls = cols[i+1].split(',')
            if len(test_calls) < min_support:
                continue
            test_alleles.append(float(cols[i]))
            supports[float(cols[i])] = len(test_calls)
        expanded_loci[locus] = test_alleles, ['-'], '-', supports

    olaps = defaultdict(list)
    for cols in common_loci:
        olaps[tuple(cols[:3])].append(cols)
    pick_best_olap(olaps)
    
    for loc in olaps.keys():
        cols = olaps[loc][0]
        test_cols = cols[:9]
        locus = (test_cols[0], int(test_cols[1]), int(test_cols[2]), test_cols[3], float(test_cols[4]))
        ref_size = locus[-1]

        expanded_alleles = []
        supports = {}
        control_alleles = []
        pvals = []
        # loop thru test allele
        for i in range(6, 9, 2):
            if test_cols[i] == '-':
                continue
            test_calls = list(map(float, test_cols[i].split(',')))
            if len(test_calls) < min_support:
                continue

            test_allele = float(test_cols[i-1])
            if float(test_cols[i-1]) - ref_size < min_expansion:
                continue

            # loop thru controls
            not_expanded = False
            control_cols = cols[9:9+10]
            control_pvals = []

            all_controls = {}
            if not from_catalog:
                for j in range(-3, 0, 2):
                    if control_cols[j] == '-':
                        continue
                    control_allele = float(control_cols[j-1])
                    if not control_allele in all_controls:
                        all_controls[control_allele] = list(map(float, control_cols[j].split(',')))
            else:
                a = 5 if use_size else 4
                alleles = [float(a) for a in control_cols[a].split(';') if a!= '-']
                #print('yy', locus, sorted(alleles), get_larger_alleles(alleles, min_count=min_support))
                alleles = get_larger_alleles([float(a) for a in control_cols[a].split(';') if a!= '-'], min_count=min_support)
                allele = np.max(alleles)
                all_controls[allele] = alleles

            control_alleles = sorted(all_controls.keys())
            for control_allele in control_alleles:
                control_calls = all_controls[control_allele]
                result = stats.ttest_ind(np.array(control_calls), np.array(test_calls), alternative='less')
                #print('uu', locus, test_calls, sorted(control_calls), result, result.pvalue)
                if result.pvalue > pval_cutoff or test_allele - control_allele + 1 < min_expansion:
                    not_expanded = True
                else:
                    control_pvals.append(np.format_float_scientific(result.pvalue, precision=2))

            if not not_expanded:
                expanded_alleles.append(test_allele)
                supports[test_allele] = len(test_calls)
                pvals.append(','.join(control_pvals))
            
        if expanded_alleles:
            expanded_loci[locus] = expanded_alleles, [','.join(list(map(str, control_alleles)))], ';'.join(pvals), supports

    return expanded_loci

def vs_all_controls(vs_controls):
    expanded_loci = {}
    for locus in set.intersection(*map(set, vs_controls)):
        expanded_alleles = set.intersection(*map(set, [vc[locus][0] for vc in vs_controls]))
        if not expanded_alleles:
            continue

        supports = []
        for allele in expanded_alleles:
            for vc in vs_controls:
                if allele in vc[locus][3]:
                    supports.append(vc[locus][3][allele])
                    break

        control_alleles = []
        pvals = []
        control_alleles_list = []
        for vc in vs_controls:
            if vc[locus][1]:
                control_alleles.append(vc[locus][1][0])
                for allele in vc[locus][1][0].split(','):
                    if allele != '-':
                        control_alleles_list.append(float(allele))
                pvals.append(vc[locus][2])
            else:
                control_alleles.append('-')
                pvals.append('-')
       
        if not control_alleles_list:
            expansion = np.median(list(expanded_alleles))
        else:
            expansion = round(np.median(list(expanded_alleles)) - np.median(control_alleles_list), 1)
        
        expanded_loci[locus] = [list(expanded_alleles), control_alleles, pvals, ','.join(map(str, list(supports))), expansion]

    return expanded_loci

def olap_gtf(bed, gtf_file):
    gtf = BedTool(gtf_file)

    skipped_types = ('retained_intron', 'nonsense_mediated_decay', 'processed_transcript')
    olaps = defaultdict(list)
    assigned = {}
    for cols in bed.intersect(gtf, wao=True, f=1.0):
        locus = tuple(cols[:3])

        if cols[3] == '.':
            assigned[locus] = None
            continue
        feature = cols[5]

        gff = create_interval_from_list(cols[3:-1])
        gene = gff.attrs['gene_name']
        transcript = None
        transcript_type = None
        gene_type = None
        if 'transcript_type' in gff.attrs and gff.attrs['transcript_type'] in skipped_types:
            continue
        if 'gene_type' in gff.attrs and gff.attrs['gene_type'] in skipped_types:
            continue

        if 'transcript_id' in gff.attrs:
            transcript = gff.attrs['transcript_id']
        if 'transcript_type' in gff.attrs:
            transcript_type = gff.attrs['transcript_type']
        if 'gene_type' in gff.attrs:
            gene_type = gff.attrs['gene_type']

        olaps[locus].append((feature, gene, gene_type))

    for locus, features in olaps.items():
        cds = [f for f in features if f[0] == 'CDS']
        exons = [f for f in features if f[0] == 'exon']
        transcripts = [f for f in features if f[0] == 'transcript']
        
        if cds:
            gene_counts = Counter(f[1] for f in cds)
            gene = gene_counts.most_common(1)[0][0]
            assigned[locus] = '{}:{}'.format('CDS', gene)
        elif exons:
            gene_counts = Counter(f[1] for f in exons)
            gene = gene_counts.most_common(1)[0][0]
            assigned[locus] = '{}:{}'.format('exon', gene)
        elif transcripts:
            coding_genes = [t for t in transcripts if t[2] == 'protein_coding']
            if coding_genes:
                gene_counts = Counter(f[1] for f in coding_genes)
                gene = gene_counts.most_common(1)[0][0]
                assigned[locus] = '{}:{}'.format('intron', gene)
            else:
                gene_counts = Counter(f[1] for f in transcripts)
                gene = gene_counts.most_common(1)[0][0]
                assigned[locus] = '{}:{}'.format('intron', gene) 

    return assigned

def olap_bed(bed, annot_file, feature_type=None):
    annot_bed = BedTool(annot_file)
    olaps = defaultdict(list)
    for cols in bed.intersect(annot_bed, wao=True, f=1.0).moveto('ee.bed'):
        locus = tuple(cols[:3])
        feature = cols[6]
        if feature == '.':
            continue
        if feature_type is not None:
            feature = '{}:{}'.format(feature_type, feature)
        olaps[locus].append(feature)

    return olaps

def locate_events(expanded_loci, gtf, enhancers_bed, promoters_bed):
    bed = create_simple_bed(list(expanded_loci.keys()))
    genic_locs = olap_gtf(bed, gtf) if gtf is not None else {}
    promoter_olaps = olap_bed(bed, promoters_bed, feature_type='promoter') if promoters_bed is not None else {}
    enhancer_olaps = olap_bed(bed, enhancers_bed, feature_type='enhancer') if enhancers_bed is not None else {}

    for locus in expanded_loci.keys():
        key = tuple(map(str, locus[:3]))
        if key in genic_locs and genic_locs[key] is not None:
            expanded_loci[locus].append(genic_locs[key])
        elif key in promoter_olaps:
            expanded_loci[locus].append(','.join(promoter_olaps[key]))
        elif key in enhancer_olaps:
            expanded_loci[locus].append(','.join(enhancer_olaps[key]))
        else:
            expanded_loci[locus].append('-')

def output(expanded_loci, out_file, has_gene=False):
    header = ['chrom', 'start', 'end', 'repeat', 'ref_size', 'test_allele', 'test_allele_support', 'expansion', 'control_alleles', 'pvals']
    if has_gene:
        header.append('gene')
    with open(out_file, 'w') as out:
        out.write('#{}\n'.format('\t'.join(header)))
        for locus in sorted(expanded_loci.keys(), key=itemgetter(0,1,2)):
            expanded_alleles = ','.join(list(map(str, expanded_loci[locus][0])))
            control_alleles = ';'.join(expanded_loci[locus][1])
            if control_alleles == '-;-':
                control_alleles = '-'
            pvals = ';'.join(expanded_loci[locus][2])

            if not re.search('\d', control_alleles):
                control_alleles = '-'
                pvals = '-'

            cols = list(map(str, locus)) + [expanded_alleles, expanded_loci[locus][3], str(expanded_loci[locus][4]), control_alleles, pvals]
            if has_gene:
                cols.append(expanded_loci[locus][5])
            out.write('{}\n'.format('\t'.join(cols)))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("test", type=str, help="genotype of test_sample")
    parser.add_argument("output", type=str, help="output")
    parser.add_argument("--controls", type=str, nargs='+', help="controls")
    parser.add_argument("--catalog", type=str, help="controls are from catalog")
    parser.add_argument("--use_copy_number", action='store_true', help="use copy number")
    parser.add_argument("--min_expansion", type=int, default=100, help="minimum expansion. Default:100")
    parser.add_argument("--min_support", type=int, default=4, help="minimum support. Default:4")
    parser.add_argument("--skip_chroms", type=str, nargs='+', help="skip chromosomes")
    parser.add_argument("--pval_cutoff", type=float, default=0.001,  help="p-value cutoff for testing T-test hypothesis. Default:0.001")
    parser.add_argument("--gtf", type=str, help="gtf")
    parser.add_argument("--promoters", type=str, help="promoters bed file")
    parser.add_argument("--enhancers", type=str, help="enhancers bed file")
    parser.add_argument("--old_version", action='store_true', help="old version of Straglr used")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    control_results = []
    if args.catalog:
        control_results.append(args.catalog)
    if args.controls:
        if len(args.controls) == 1 and os.path.exists(args.controls[0]) and os.path.splitext(args.controls[0])[1] != '.tsv':
            with open(args.controls[0], 'r') as ff:
                control_results.extend([f.rstrip() for f in ff.readlines() if os.path.exists(f.rstrip()) and os.path.splitext(f.rstrip())[1] == '.tsv'])
        else:
            control_results.extend([c for c in args.controls if os.path.exists(c)])

    if not control_results:
        sys.exit()

    expanded_loci = {}
    use_size = not args.use_copy_number
    test_bed = parse_straglr_tsv(args.test, use_size=use_size, skip_chroms=args.skip_chroms, old_version=args.old_version)

    if control_results:
        vs_controls = []
        for i in range(len(control_results)):
            control_result = control_results[i]
            print('comparing {} vs {}'.format(args.test, control_result))
            from_catalog = False
            if i==0 and args.catalog:
                from_catalog = True
                control_bed = control_result
            else:
                control_bed = parse_straglr_tsv(control_result, use_size=use_size, old_version=args.old_version)
   
            vs_controls.append(vs_each_control(test_bed,
                                               control_bed,
                                               args.pval_cutoff,
                                               min_expansion=args.min_expansion,
                                               min_support=args.min_support, 
                                               label=control_result,
                                               from_catalog=from_catalog,
                                               use_size=use_size))
    
        if vs_controls:
            expanded_loci = vs_all_controls(vs_controls)

    if expanded_loci:
        has_gene = False
        if args.gtf or args.enhancers:
            has_gene = True
            locate_events(expanded_loci, args.gtf, args.enhancers, args.promoters)
        output(expanded_loci, args.output, has_gene=has_gene)

main()
