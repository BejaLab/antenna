#!/usr/bin/env python3

from Bio import SeqIO
from sys import stderr
import csv
import re
from crc64iso.crc64iso import crc64

info_file = str(snakemake.input['info'])
a2m_file = str(snakemake.input['a2m'])
pos_file = str(snakemake.input['pos'])
usearch_file = str(snakemake.input['usearch'])
blasso_file  = str(snakemake.input['blasso'])
out_file = str(snakemake.output)

infos = dict()
data_type = 'JGI'
header = []
id_cols = { 'OM-RGC_ID', 'genome_id', 'gene_oid', 'ID' }
with open(info_file) as fh:
    tsv = csv.reader(fh, delimiter = '\t')
    header = next(tsv)
    header = [v + str(header[:i].count(v) + 1) if header.count(v) > 1 else v for i, v in enumerate(header)]
    id_cols = set(header) & id_cols
    assert len(id_cols) > 0, "Id column not found"
    id_col = id_cols.pop()
    has_sequence = 'sequence' in header
    for line in tsv:
        row = dict(zip(header, line))
        record_id = row.pop(id_col).replace(' ', '_')
        if has_sequence:
            del row['sequence']
        infos[record_id] = row
    header.remove(id_col)
    if has_sequence:
        header.remove('sequence')

positions = []
groups = []
with open(pos_file) as fh:
    for line in fh:
        pos, group = line.rstrip().split()
        pos0 = int(pos) - 1
        positions.append((pos0, group))
        if len(group) > 1 and group not in groups:
            groups.append(group)

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
    tsv.writerow([ 'record_id', 'checksum', 'target', 'family', 'ident', 'wl_mean', 'wl_sd' ] + groups + header)
    with open(a2m_file) as a2m_fh:
        a2m = SeqIO.parse(a2m_fh, 'fasta')
        ref = next(a2m)
        ref_pos = 0
        aln_pos = []
        for pos in range(len(ref.seq)):
            if ref.seq[pos] not in [ '-', '.' ]:
                for pos0, group in positions:
                    if ref_pos == pos0:
                        aln_pos.append((pos, group))
                ref_pos += 1
        for rec in a2m:
            if rec.id in matches:
                seq = re.sub('[^A-Z]', '', str(rec.seq))
                cks = crc64(seq)
                match  = matches[rec.id]
                blasso = blassos[rec.id]
                info_id, *rest = rec.id.split('/')
                info   = infos[info_id]
                res = { group: [] for group in groups }
                for pos, group in aln_pos:
                    if group in groups:
                        res[group].append(rec.seq[pos])
                    elif rec.seq[pos] != group:
                        break
                else:
                    row = [ rec.id, cks, match['target'], match['family'], match['ident'], blasso['wl_mean'], blasso['wl_sd'] ]
                    row += [ ''.join(res[x]) for x in groups ]
                    row += info.values()
                    tsv.writerow(row)
