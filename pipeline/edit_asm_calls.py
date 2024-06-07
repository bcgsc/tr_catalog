from pybedtools import BedTool
import argparse
import pysam
import os
import sys
from collections import defaultdict
from operator import itemgetter
import pandas as pd
import re
import tempfile

def find_olaps(straglr_bed, asm_beds, st_ncols=10, asm_ncols=9, f=0.8):
    olaps = {}
    haps = set()
    for cols in straglr_bed.intersect(asm_beds, f=f, r=True, wao=True):
        st_call = tuple(cols[:st_ncols])
        asm_call = tuple(cols[st_ncols + 1 : st_ncols + asm_ncols])
        if asm_call[0] == '.':
            continue
        
        # hap = paternal or maternal
        hap = cols[-2].split('.')[1]
        haps.add(hap)

        if not st_call in olaps:
            olaps[st_call] = defaultdict(list)
        olaps[st_call][hap].append(asm_call)

    return olaps, find_non_olaps(straglr_bed, asm_beds, st_ncols), haps

def find_non_olaps(straglr_bed, asm_beds, ncols):
    non_olaps = []
    for cols in straglr_bed.intersect(asm_beds, v=True):
        non_olaps.append(tuple(cols[:ncols]))
    return non_olaps

def overlap_uniques(bed1, bed2, ncols=10, f=0.8):
    olaps = []
    for cols in bed1.intersect(bed2, f=f, wao=True):
        st_call1 = tuple(cols[:ncols])
        st_call2 = tuple(cols[ncols+1:ncols+1 + ncols])
        if st_call2[0] == '.':
            continue
        olaps.append((st_call1, st_call2))

    for cols in bed1.intersect(bed2, F=f, wao=True):
        st_call1 = tuple(cols[:ncols])
        st_call2 = tuple(cols[ncols+1:ncols+1 + ncols])
        if st_call2[0] == '.' or (st_call1, st_call2) in olaps:
            continue
        olaps.append((st_call1, st_call2))
 
    return olaps

def find_novels(uniques, skip_chroms=None):
    beds = []
    read_types = []
    for read_type, calls in uniques:
        if len(calls) == 0:
            continue
        bed_str = ''
        for call in calls:
            if call[0] in skip_chroms:
                continue
            bed_str += '{}\n'.format('\t'.join(list(call) + [read_type]))
        bed = BedTool(bed_str, from_string=True)
        beds.append(BedTool(bed_str, from_string=True))
        read_types.append(read_type)

    novel_calls = defaultdict(list)
    log = ''
    olaps = ([], [])
    # 2 read types
    if len(uniques) == 2:
        if 'np' in read_types:
            np_index = read_types.index('np')
        if 'hifi' in read_types:
            hifi_index = read_types.index('hifi')
        if len(beds) == 2:
            olaps = overlap_uniques(beds[0], beds[1])
            for calls in olaps:
                alleles = [extract_st_alleles(call) for call in calls]

                passed_supports = [screen_novel_call(alleles[i], read_types[i]) for i in range(len(alleles))]
                if passed_supports[hifi_index]:
                    log += 'novel\t{}\t"{}"\t"{}"\n'.format(read_type, ' '.join(calls[hifi_index]), ';'.join([','.join(list(map(str, a))) for a in alleles]))
                    novel_calls[read_types[hifi_index]].append((calls[hifi_index], alleles[hifi_index]))
                elif passed_supports[np_index]:
                    log += 'novel\t{}\t"{}"\t"{}"\n'.format(read_type, ' '.join(calls[np_index]), ';'.join([','.join(list(map(str, a))) for a in alleles]))
                    novel_calls[read_types[np_index]].append((calls[np_index], alleles[np_index]))

    for read_type, calls in uniques:
        if not calls:
            continue
        index = read_types.index(read_type)
        used = []
        if olaps[0]:
            used = [o[index] for o in olaps]
        for call in calls:
            if not call in used and not call[0] in skip_chroms:
                allele = extract_st_alleles(call)
                #print('rr', read_type, call, screen_novel_call(allele, read_type))
                if screen_novel_call(allele, read_type):
                    log += 'novel\t{}\t"{}"\t"{}"\n'.format(read_type, ' '.join(call), ';'.join([','.join(list(map(str, a))) for a in allele]))
                    novel_calls[read_type].append((call, allele))

    novel_entries = []
    for read_type in novel_calls.keys():
        for call, alleles in novel_calls[read_type]:
            for allele, i in zip(alleles, (0,1)):
                size_change = int(allele[0]) - (int(call[2]) - int(call[1]) + 1)
                cols = list(call[:4]) + ['straglr_{}'.format(read_type)] + list(allele[:2]) + [size_change]
                novel_entries.append((i, (), tuple(map(str, cols))))

    totals = {}
    for read_type in novel_calls.keys():
        totals['novel_{}'.format(read_type)] = len(novel_calls[read_type])
    
    return novel_entries, log, totals

def screen_novel_call(alleles, read_type):
    if len(alleles) == 1:
        min_support = {(0,5000):{'np':10, 'hifi':10}, (5000,1000000):{'np':4, 'hifi':4}}
    elif len(alleles) == 2:
        min_support = {(0,5000):{'np':10, 'hifi':10}, (5000,1000000):{'np':4, 'hifi':4}}

    passed = [a for a in alleles if screen_by_support(a[0], a[2], min_support, read_type)]
    return len(alleles) == len(passed)
        
def extract_st_alleles(cols):
    alleles = []
    for i in range(4, len(cols), 3):
        # size copy_number support
        if cols[i] != '-':
            allele = tuple(map(float, cols[i:i+3]))
            alleles.append(allele)

    # remove alleles with few supporting reads
    alleles_sorted = sorted(alleles, key=itemgetter(0), reverse=True)
    if len(alleles_sorted) == 2 and alleles_sorted[-1][2] < max(4, 0.1 * alleles_sorted[0][2]):
        alleles = [alleles_sorted[0]]

    # homozygous
    if len(alleles) == 1:
        alleles.append(alleles[0])

    return sorted(alleles, key=itemgetter(0), reverse=True)

def agree(a1, a2, f=0.1):
    d = abs(a1-a2)
    return  d <= f * a1 or d <= f * a2

def find_matches(calls1, calls2):
    matched = []
    used1 = []
    used2 = []
    for i in range(len(calls1)):
        if i in used1:
            continue
        for j in range(len(calls2)):
            if i in used1:
                break
            if j in used2:
                continue
            if agree(calls1[i], calls2[j]):
                matched.append(calls1[i])
                used1.append(i)
                used2.append(j)
    return matched, [calls2[j] for j in used2]

def screen_by_support(size, support, min_support, tech):
    for size_range in min_support:
        if size >= size_range[0] and size <= size_range[1]:
            return support >= min_support[size_range][tech]
    return False

def check_asm_calls(np_olaps, hifi_olaps, haps, min_size=1000):
    validated = defaultdict(list)
    unvalidated = defaultdict(list)
    all_matches = defaultdict(dict)
    all_asm_calls = set()
    read_types = ('np', 'hifi')

    totals = {}
    for read_type, olaps in zip(read_types, (np_olaps, hifi_olaps)):
        if not read_type in totals:
            totals['{}_{}'.format(read_type, 'olapped')] = 0
            totals['{}_{}'.format(read_type, 'all_matched')] = 0 
            totals['{}_{}'.format(read_type, 'unmatched')] = 0
        for st_call in olaps.keys():
            st_alleles = extract_st_alleles(st_call)
            st_sizes = [a[0] for a in st_alleles]

            asm_call = None
            if haps[0] in olaps[st_call] and haps[1] in olaps[st_call]:
                if len(olaps[st_call][haps[0]]) == 1 and len(olaps[st_call][haps[1]]) == 1:
                    asm_call = olaps[st_call][haps[0]][0], olaps[st_call][haps[1]][0]
            elif haps[0] in olaps[st_call] and len(olaps[st_call][haps[0]]) == 1:
                asm_call = (olaps[st_call][haps[0]][0], ())
            elif haps[1] in olaps[st_call] and len(olaps[st_call][haps[1]]) == 1:
                asm_call = ((), olaps[st_call][haps[1]][0])
            
            if asm_call is not None:
                totals['{}_{}'.format(read_type, 'olapped')] += 1
                asm_sizes = [float(a[5]) for a in asm_call if a]

                asm_sizes_matched, st_sizes_matched = find_matches(asm_sizes, st_sizes)
                asm_sizes_unmatched = set(asm_sizes) - set(asm_sizes_matched)
                st_sizes_unmatched = set(st_sizes) - set(st_sizes_matched)

                all_matches[asm_call][read_type] = (st_call, asm_sizes, st_sizes, asm_sizes_matched, st_sizes_matched, asm_sizes_unmatched, st_sizes_unmatched)
                if len(asm_sizes_matched) == len(asm_call):
                    validated[asm_call].append(read_type)
                    totals['{}_{}'.format(read_type, 'all_matched')] += 1
                else:
                    #print('tt', read_type, st_call, st_sizes, asm_call, asm_sizes, st_sizes_matched, asm_sizes_matched, st_sizes_unmatched, asm_sizes_unmatched)
                    unvalidated[asm_call].append(read_type)
                    totals['{}_{}'.format(read_type, 'unmatched')] += 1

    corrections = []
    log = ''
    num_edits = 0
    for asm_call in unvalidated.keys():
        # unvalidated by both read types
        if len(unvalidated[asm_call]) == 2:
            results = check_both(asm_call, all_matches[asm_call], haps)
            corrections.extend(results[0])
            log += results[1]
            if results[0]:
                num_edits += 1
        else:
            read_type = unvalidated[asm_call][0]
            results = check_single(asm_call, read_type, all_matches[asm_call][read_type], haps)
            corrections.extend(results[0])
            log += results[1]
            if results[0]:
                num_edits += 1
    
    totals['edits'] = num_edits
    return corrections, log, totals

def check_single(asm_call, read_type, matches, haps):
    min_support = {(0,5000):{'np':10, 'hifi':10}, (5000,1000000):{'np':4, 'hifi':4}}
    # for nanopore
    min_size = {'hifi':0, 'np':1000}

    st_call = matches[0]
    st_alleles = extract_st_alleles(st_call)
    st_support = dict((float(st_call[i]), int(st_call[i+2])) for i in (4,7) if st_call[i] != '-')
    
    asm_sizes_matched, st_sizes_matched, asm_sizes_unmatched, st_sizes_unmatched = matches[-4:]
    asm_sizes = []
    for a in asm_call:
        asm_sizes.append(float(a[5]) if len(a) > 0 else None)

    #print('ee', asm_call, st_call, read_type, st_support, asm_sizes_matched, st_sizes_matched, asm_sizes_unmatched, st_sizes_unmatched, asm_sizes)

    corrections = []
    log = ''
    if asm_sizes_matched and st_sizes_unmatched:
        if asm_sizes_unmatched:
            # don't check valdiated because hifi may validate both asm alleles but there is a valid np allele
            if len(st_sizes_unmatched) == 1 and len(asm_sizes_unmatched) == 1:
                st_size_unmatched = list(st_sizes_unmatched)[0]
                asm_size_unmatched = list(asm_sizes_unmatched)[0]
                ms = min_size[read_type] if read_type == 'hifi' else max(min_size[read_type], asm_size_unmatched)

                if not agree(st_size_unmatched, asm_size_unmatched):
                    # change regardless of whether np confirms or not because np more messy? or only if it is not confirmed by other?
                    # should be former, logic for changing muc3a allele
                    if screen_by_support(st_size_unmatched, st_support[st_size_unmatched], min_support, read_type) and st_size_unmatched >= ms:
                        #new_size = st_size_unmatched
                        unmatched_st_allele = [s for s in st_alleles if s[0] == st_size_unmatched][0]
                        ref_size = int(st_call[2]) - int(st_call[1]) + 1
                        new_allele = list(st_call[:4]) + ['straglr_{}'.format(read_type)] + list(unmatched_st_allele[:2]) + [ref_size]

                        i = asm_sizes.index(asm_size_unmatched)
                        #hap = 'paternal' if i == 0 else 'maternal'
                        hap = haps[0] if i == 0 else haps[1]
                        #print('ww1', asm_call, asm_sizes, read_type, asm_sizes_matched, st_sizes_matched, st_call, i, asm_call[i], hap, new_allele, unmatched_st_allele)
                        log += 'correct\t{}\t{}\t"{}"\t"{}"\t"{}"\n'.format(read_type, hap, ' '.join(asm_call[i]), ' '.join(st_call), ' '.join(list(map(str, unmatched_st_allele))))
                        corrections.append((i, asm_call[i], new_allele))

        # homozygous asm and both matched, but not straglr (e.g. MUC3A)
        else:
            st_size_unmatched = list(st_sizes_unmatched)[0]
            asm_size_matched = list(asm_sizes_matched)[0]
            ms = min_size[read_type] if read_type == 'hifi' else max(min_size[read_type], asm_size_matched)
            if screen_by_support(st_size_unmatched, st_support[st_size_unmatched], min_support, read_type):
            #if screen_by_support(st_size_unmatched, st_support[st_size_unmatched], min_support, read_type) and st_size_unmatched >= ms:
                #new_size = st_size_unmatched
                unmatched_st_allele = [s for s in st_alleles if s[0] == st_size_unmatched][0]
                ref_size = int(st_call[2]) - int(st_call[1]) + 1
                new_allele = list(st_call[:4]) + ['straglr_{}'.format(read_type)] + list(unmatched_st_allele[:2]) + [ref_size]
                
                # arbitrary pick paternal allele
                i = 0 if asm_sizes.index(asm_size_matched) == 1 else 1
                hap = haps[0] if i == 0 else haps[1]
                #print('ww2', asm_call, asm_sizes, read_type, asm_sizes_matched, st_sizes_matched, st_call, i, asm_call[i], hap, new_allele, unmatched_st_allele)
                log += 'correct\t{}\t{}\t{}\t"{}"\t"{}"\n'.format(read_type, hap, "missing", ' '.join(st_call), ' '.join(list(map(str, unmatched_st_allele))))
                corrections.append((i, asm_call[i], new_allele))

    return corrections, log

def check_both(asm_call, matches, haps):
    corrections = []
    log = ''
    asm_calls_edited = set()
    st_calls_used = set()
    read_types = list(matches.keys())
    # same asm call unvalidated by both np and hifi
    if matches[read_types[0]][-2] == matches[read_types[1]][-2]:
        #print('dd', asm_call, matches[read_types[0]][0], matches[read_types[1]][0], matches[read_types[0]][-2], matches[read_types[1]][-2])

        # same number of unmatched straglr calls for both np and hifi (>0)
        if matches[read_types[0]][-1] and len(matches[read_types[0]][-1]) == len(matches[read_types[1]][-1]):
            unmatched_asm_size = None
            if matches[read_types[0]][-2]:
                unmatched_asm_size = list(matches[read_types[0]][-2])[0]

            if len(matches[read_types[0]][-1]) == 1:
                st_sizes = [list(matches[r][-1])[0] for r in read_types]
                #print('xx', asm_call, len(matches[read_types[0]][-1]), len(matches[read_types[0]][-1]), unmatched_asm_size, st_sizes)
                # same unmatched st allele from both np and hifi
                if agree(st_sizes[0], st_sizes[1]):
                    #new_size = st_sizes[1]
                    # used hifi arbitarily
                    st_call = matches[read_types[1]][0]
                    st_alleles = extract_st_alleles(st_call)
                    unmatched_st_allele = [s for s in st_alleles if s[0] == st_sizes[1]][0]
                    ref_size = int(st_call[2]) - int(st_call[1]) + 1
                    new_allele = list(st_call[:4]) + ['straglr_hifi'] + list(unmatched_st_allele[:2]) + [ref_size]

                    for i in range(len(asm_call)):
                        # either missing one allele in asm, or homozygous(both matched) just arbitrarily pick paternal allele or find the unmatched allele
                        if len(asm_call[i]) == 0 or (unmatched_asm_size is None and i==0) or (unmatched_asm_size is not None and float(asm_call[i][5]) == unmatched_asm_size):
                            hap = haps[0] if i == 0 else haps[1]
                            ac = '"{}"'.format(' '.join(asm_call[i])) if len(asm_call[i]) > 0 else 'missing'
                            log += 'correct\t{}\t{}\t{}\t"{}"\t"{}"\n'.format(read_types[1], hap, ac, ' '.join(st_call), ' '.join(list(map(str, unmatched_st_allele)))) 
                            #print('bb1', asm_call, i, hap, new_allele)
                            corrections.append((i, asm_call[i], new_allele))
                            asm_calls_edited.add(tuple(asm_call[i]))
                            st_calls_used.add(tuple(st_call))
            else:
                # 2 unmatched hifi and np alleles (implies 0 matched asm alleles
                st_sizes1 = sorted(list(matches[read_types[0]][-1]))
                st_sizes2 = sorted(list(matches[read_types[1]][-1]))
                #print('kk', asm_call, matches[read_types[0]][0], matches[read_types[1]][0], matches[read_types[0]][-1], matches[read_types[1]][-1], st_sizes1, st_sizes2, agree(st_sizes1[0], st_sizes2[0]), agree(st_sizes1[1], st_sizes2[1]))
                if agree(st_sizes1[0], st_sizes2[0]) and agree(st_sizes1[1], st_sizes2[1]):
                    #print('kk1', asm_call, matches[read_types[0]][0], matches[read_types[1]][0], matches[read_types[0]][-1], matches[read_types[1]][-1])
                    st_call = matches[read_types[1]][0]
                    st_alleles = extract_st_alleles(st_call)
                    for i in range(len(asm_call)):
                        ref_size = int(st_call[2]) - int(st_call[1]) + 1
                        new_allele = list(st_call[:4]) + ['straglr_hifi'] + list(st_alleles[i][:2]) + [ref_size]
                        hap = haps[0] if i == 0 else haps[1]
                        ac = '"{}"'.format(' '.join(asm_call[i])) if len(asm_call[i]) > 0 else 'missing'
                        log += 'correct\t{}\t{}\t{}\t"{}"\t"{}"\n'.format(read_types[1], hap, ac, ' '.join(st_call), ' '.join(list(map(str, st_alleles[i]))))
                        corrections.append((i, asm_call[i], new_allele))
                        asm_calls_edited.add(tuple(asm_call[i]))
                        st_calls_used.add(tuple(st_call))
                #else:
                    #print('kk2', asm_call, matches[read_types[0]][0], matches[read_types[1]][0], matches[read_types[0]][-1], matches[read_types[1]][-1])

    return corrections, log, {'asm_calls_edited': len(asm_calls_edited), 'st_calls_used': len(st_calls_used)}

def edit_bed(in_bed, out_bed, remove=[], add=[]):
    header = ['chrom', 'start', 'end', 'motif', 'seq', 'size', 'copy_number', 'size_change', 'hap']
    fd, tmp_tsv = tempfile.mkstemp()
    
    BedTool(in_bed).saveas(tmp_tsv)
    df = pd.read_csv(tmp_tsv, sep='\t', names=header)
    hap = df['hap'].iloc[0]

    index_remove = []
    for r in remove:
        index_remove = df[(df.chrom == r[0]) & (df.start == int(r[1])) & (df.end == int(r[2]))].index
        df.drop(index_remove, inplace=True)
    
    df_new = pd.DataFrame(columns=header)
    for i in range(len(add)):
        df_new.loc[i] = list(add[i]) + [hap]

    df_out = pd.concat([df, df_new])
    df_out.to_csv(out_bed, index=False, sep='\t', header=False)
    BedTool(out_bed).tabix(in_place=True, force=True)

    os.remove(tmp_tsv)

def make_changes(corrections, in_beds, out_beds):
    for i in range(len(in_beds)):
        replace = [(c[1], c[2]) for c in corrections if c[0] == i and len(c[1]) > 0]
        add = [c[2] for c in corrections if c[0] == i and len(c[1]) == 0]

        remove = [tuple(r[0]) for r in replace]
        add.extend([tuple(r[1]) for r in replace])

        bed = edit_bed(in_beds[i], out_beds[i], remove=remove, add=add)

def output_log(out_file, log, totals):
    with open(out_file, 'w') as out:
        for cat, num in totals.items():
            out.write('#{}: {}\n'.format(cat, num))
        out.write(log)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("asms", nargs=2, type=str, help="haploid asm results")
    parser.add_argument("out_suffix", type=str, help="output suffix")
    parser.add_argument("--np", type=str, help="nanopore straglr bed")
    parser.add_argument("--hifi", type=str, help="hifi straglr bed")
    parser.add_argument("--log", type=str, help="log file")
    parser.add_argument("--skip_chroms", type=str, default="chrY", help="comma-separated list of skip chromosomes, e.g. chrY")
    parser.add_argument("--cen_bed", type=str, help="centromere bed")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    skip_chroms = args.skip_chroms.split(',')

    np_olaps, np_unique = {}, {}
    hifi_olaps, hifi_unique = {}, {}
    haps = []
    if args.np:
        np_bed = BedTool(args.np)
        if args.cen_bed:
            np_bed = np_bed.intersect(args.cen_bed, v=True)
        np_olaps, np_unique, haps = find_olaps(np_bed, args.asms)
    if args.hifi:
        hifi_bed = BedTool(args.hifi)
        if args.cen_bed:
            hifi_bed = hifi_bed.intersect(args.cen_bed, v=True)
        hifi_olaps, hifi_unique, haps = find_olaps(hifi_bed, args.asms)
  
    if len(haps) != 2:
        sys.exit('not getting 2 haps:', haps)
    elif 'paternal' in haps and 'maternal' in haps:
        haps = ('paternal', 'maternal')
    else:
        haps = sorted(list(haps))
    print('haps', haps)

    corrections, log_corrections, totals_edits = check_asm_calls(np_olaps, hifi_olaps, haps)
    additions, log_additions, totals_novels = find_novels((('np', np_unique), ('hifi', hifi_unique)), skip_chroms=skip_chroms)

    totals = totals_edits
    totals.update(totals_novels)

    changes = corrections + additions
    if changes:
        in_beds = args.asms if haps[0] in args.asms[0] else args.asms[::-1]
        out_beds = []
        for fn in in_beds:
            fn_prefix, fn_ext = os.path.splitext(re.sub('.gz$', '', fn))
            out_beds.append('{}.{}{}'.format(fn_prefix, args.out_suffix, fn_ext))
        make_changes(changes, in_beds, out_beds)
        
        output_log(args.log, log_corrections + log_additions, totals)

main()
