#!/usr/bin/env python3

from Bio import SeqIO
from sys import stderr
import csv
import re

info_file = str(snakemake.input['info'])
a2m_file = str(snakemake.input['a2m'])
pos_file = str(snakemake.input['pos'])
usearch_file = str(snakemake.input['usearch'])
blasso_file  = str(snakemake.input['blasso'])
out_file = str(snakemake.output)

infos = dict()
data_type = 'JGI'
header = []
with open(info_file) as fh:
    tsv = csv.reader(fh, delimiter = '\t')
    header = next(tsv)
    header = [v + str(header[:i].count(v) + 1) if header.count(v) > 1 else v for i, v in enumerate(header)]
    if 'OM-RGC_ID' in header:
        data_type = 'OM-RGC'
    elif 'genome_id' in header:
        data_type = 'GEM'
    for line in tsv:
        row = dict(zip(header, line))
        if data_type == 'OM-RGC':
            record_id = row.pop('OM-RGC_ID')
            del row['sequence']
        elif data_type == 'GEM':
            record_id = row.pop('genome_id')
        else:
            record_id = row.pop('gene_oid').replace(' ', '_')
        infos[record_id] = row
    if data_type == 'OM-RGC':
        header.remove('OM-RGC_ID')
        header.remove('sequence')
    elif data_type == 'GEM':
        header.remove('genome_id')
    else:
        header.remove('gene_oid')

positions = dict()
with open(pos_file) as fh:
    for line in fh:
        pos, *expected = line.rstrip().split()
        pos0 = int(pos) - 1
        if expected:
            positions[pos0] = expected[0]
        else:
            positions[pos0] = ''

blassos = dict()
with open(blasso_file) as fh:
    for line in fh:
        record_id, wl_mean, wl_sd = line.split()
        blassos[record_id] = { 'wl_mean': wl_mean, 'wl_sd': wl_sd }

matches = dict()
with open(usearch_file) as fh:
    for line in fh:
        query, target, ident = line.split('\t')
        if query not in matches:
            target_id, *rest = target.split()
            search = re.search('{(.+?)}', target)
            assert search, "Unexpected reference header format: %s" % target
            family = search.group(1)
            matches[query] = { 'target': target, 'family': family, 'ident': float(ident) }

aminoacids = [ 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' ]
with open(out_file, 'w') as out_fh:
    tsv = csv.writer(out_fh, delimiter = "\t")
    tsv.writerow([ 'record_id', 'positions', 'target', 'family', 'ident', 'wl_mean', 'wl_sd' ] + header)
    with open(a2m_file) as a2m_fh:
        a2m = SeqIO.parse(a2m_fh, 'fasta')
        ref = next(a2m)
        ref_pos = 0
        aln_pos = dict()
        for pos in range(len(ref.seq)):
            if ref.seq[pos] not in [ '-', '.' ]:
                if ref_pos in positions:
                    assert ref.seq[pos] in aminoacids, "Reference position %d is not aligned" % (ref_pos + 1)
                    aln_pos[pos] = positions[ref_pos]
                ref_pos += 1
        for rec in a2m:
            if rec.id in matches:
                match  = matches[rec.id]
                blasso = blassos[rec.id]
                info_id, *rest = rec.id.split('/')
                info   = infos[info_id]
                res = []
                for pos, expected in aln_pos.items():
                    if rec.seq[pos] not in aminoacids:
                        break
                    elif not expected:
                        res.append(rec.seq[pos])
                    elif rec.seq[pos] != expected:
                        break
                else:
                    row = [ rec.id, ''.join(res), match['target'], match['family'], match['ident'], blasso['wl_mean'], blasso['wl_sd'] ]
                    row += info.values()
                    tsv.writerow(row)