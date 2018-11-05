import sys
import os
import gzip
import math
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../align/')
import align_helpers
sys.path.insert(0, '../misc/')
import read_fasta
import read_table

def calculate_polymorphism(seqs):
    poly_sites = []
    keys = seqs.keys()
    n = len(seqs[seqs.keys()[0]])
    ignore = set([gp.unsequenced_symbol.upper(), gp.unsequenced_symbol.lower()])
    t = 0
    for i in range(n):
        a = [seqs[key][i] for key in keys]
        if gp.gap_symbol in a:
            continue
        a_set = set(a) - ignore
        if len(a_set) > 1:
            poly_sites.append(i)
        if len(a_set) > 0:
            t += 1
    return len(poly_sites), t

def diffs_per_site(s1, s2):
    ignore = set([gp.unsequenced_symbol.lower(), gp.gap_symbol])
    num = 0
    den = 0
    for i in range(len(s1)):
        if s1[i] in ignore or s2[i] in ignore:
            continue
        den += 1
        if s1[i] != s2[i]:
            num += 1
    if den == 0:
        return 'NA'
    return float(num)/den

def calculate_nuc_div(seqs):
    # average number of nucleotide differences per site between a pair
    # of sequences
    keys = seqs.keys()
    n = len(keys)
    num = 0
    den = 0
    for i in range(n-1):
        for j in range(i, n):
            d = diffs_per_site(seqs[keys[i]].lower(), seqs[keys[j]].lower())
            if d != 'NA':
                num += d
                den += 1
    if den == 0:
        return 'NA'
    return float(num)/den

# read in shared regions
shared_regions , l = \
    read_table.read_table_rows('shared_introgression_nonsingleton_list.txt', '\t')

# read in strain dirs information
s = align_helpers.get_strains(align_helpers.flatten(gp.non_ref_dirs.values()))
strain_dirs = dict(s)

# for each shared region:
# - calculate fraction of sites that are polymorphic among introgressed strains
# - for each introgressed strain, calculate:
#   - number of unique variants among introgressed strains (or all strains?)
f = open('shared_introgression_nonsingleton_polymorphism.txt', 'w')
f.write('region_number\tchromosome\tstart\tend\tpi\tfrac_poly\tnum_poly\tnum_total\tnum_strains\tstrain_list\n')
for chrm in gp.chrms:

    chrom_seqs = {}
    for region_number in shared_regions.keys():
        if shared_regions[region_number]['chromosome'] != chrm:
            continue
        
        print '*', region_number

        strains = shared_regions[region_number]['strain_list'].split(',')
        ref_start = int(shared_regions[region_number]['start'])
        ref_end = int(shared_regions[region_number]['end'])
    
        region_seqs = {}
        for strain in strains:
            print ' ', strain
            ref_to_strain_coords = [float(x[:-1]) for x in \
                                    gzip.open(gp.analysis_out_dir_absolute + \
                                              'coordinates/S288c_to_' + strain + \
                                              '_chr' + chrm + '.txt.gz').readlines()]
            strain_start = int(max(0, math.ceil(ref_to_strain_coords[ref_start])))
            strain_end = int(math.floor(ref_to_strain_coords[ref_end]))
            
            if not chrom_seqs.has_key(strain):
                chrom_seqs[strain] = read_fasta.read_fasta(strain_dirs[strain] + \
                                                           strain + '_chr' + \
                                                           chrm + gp.fasta_suffix)[1][0]
            #seq = chrom_seqs[strain][strain_start:strain_end+1]
            seq = [gp.gap_symbol for i in range(ref_start, ref_end + 1)]
            for i in range(ref_start, ref_end + 1):
                c = ref_to_strain_coords[i]
                if int(c) == c:
                    seq[i - ref_start] = chrom_seqs[strain][int(c)]
            region_seqs[strain] = ''.join(seq)

        p, t = calculate_polymorphism(region_seqs)
        fp = 'NA'
        if t != 0:
            fp = float(p)/t

        nuc_div = calculate_nuc_div(region_seqs)
        
        f.write(region_number + '\t')
        f.write(chrm + '\t')
        f.write(shared_regions[region_number]['start'] + '\t')
        f.write(shared_regions[region_number]['end'] + '\t')
        f.write(str(nuc_div) + '\t')
        f.write(str(fp) + '\t')
        f.write(str(p) + '\t')
        f.write(str(t) + '\t')
        f.write(shared_regions[region_number]['num_strains'] + '\t')
        f.write(shared_regions[region_number]['strain_list'] + '\n')
        f.flush()
f.close()
