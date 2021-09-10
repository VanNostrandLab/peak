#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Utility for plotting peaks in a certain region.
"""

import os
import sys
import math
import argparse
import itertools

import cmder
import inflect
import pandas as pd
import pyBigWig
from seqflow import Flow, task, logger
import matplotlib.pylab as plt
import seaborn as sns
from matplotlib import rcParams
import matplotlib.colors

rcParams['font.sans-serif'] = ['Times New Roman', 'Tahoma', 'DejaVu Sans', 'Lucida Grande', 'Verdana']
plt.switch_backend('agg')

parser = argparse.ArgumentParser(description=__doc__, prog='peak-plot')
parser.add_argument('--peak_beds', nargs='+', help="Space separated peak bed files.")
parser.add_argument('--pos_bigwigs', nargs='+', help="Path to BigWig files contains data on positive strand.")
parser.add_argument('--neg_bigwigs', nargs='+', help="Path to BigWig files contains data on negative strand.")
parser.add_argument('--outdir', type=str, help="Path to output directory, default: current work directory.")
parser.add_argument('--bases', type=int, help="Number of bases need to expand from peak center, default: 50.",
                    default=50)
parser.add_argument('--peak_bed_ids', nargs='+', help="Optional space separated short IDs (e.g., S1, S2, ...) for "
                                                      "peak bed files, default: basename of each peak bed file "
                                                      "without .bed extension.")
parser.add_argument('--bigwig_ids', nargs='+', help="Optional space separated short IDs (e.g., S1, S2, ...) for "
                                                    "BigWig files, default: basename of each BigWig file "
                                                    "without .bw extension.")
parser.add_argument('--cleanup', action='store_true', help='Clean up temporary files and links after complete.')
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and inputs/outputs of each step without actually running the pipeline.')
parser.add_argument('--debug', action='store_true', help='Invoke debug mode (only for develop purpose).')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

NEED_TO_DELETE = []


def right_replace(s, src, tar):
    if s.endswith(src):
        return f'{s[:-len(src)]}{tar}'
    return s


def validate_paths():
    def files_exist(files, tag):
        if not files:
            logger.error(f'No {tag} were provided, aborted.')
            sys.exit(1)
        engine, paths = inflect.engine(), []
        for i, file in enumerate(files, start=1):
            if os.path.exists(file):
                if not os.path.isfile(file):
                    logger.error(f'The {engine.ordinal(i)} file in {tag} "{file}" is not a file.')
                    sys.exit(1)
                else:
                    paths.append(os.path.abspath(file))
            else:
                logger.error(f'The {engine.ordinal(i)} file in {tag} "{file}" does not exist.')
                sys.exit(1)
        return paths

    def link_file(file, link):
        if not os.path.exists(link):
            os.symlink(file, link)
            NEED_TO_DELETE.append(link)
        return link

    peak_beds = files_exist(args.peak_beds, 'Peak beds')
    pos_bigwigs = files_exist(args.pos_bigwigs, 'BigWig positive')
    neg_bigwigs = files_exist(args.neg_bigwigs, 'BigWig negative')
    outdir = args.outdir or os.getcwd()
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            logger.error(f'Outdir "{outdir}" is a file not a directory.')
            sys.exit(1)
    else:
        logger.error(f'Outdir "{outdir}" does not exist, try to create ...')
        os.mkdir(outdir)
        logger.error(f'Successfully created Outdir "{outdir}".')

    beds, positive_bigwigs, negative_bigwigs = {}, {}, {}
    peak_bed_ids = args.peak_bed_ids if args.peak_bed_ids else [''] * len(peak_beds)
    for i, (peak_bed, basename) in enumerate(zip(peak_beds, peak_bed_ids), start=1):
        basename = basename if basename else right_replace(os.path.basename(peak_bed), '.bed', '')
        peak_bed_link = link_file(peak_bed, os.path.join(outdir, f'{basename}.peak.bed'))
        beds[basename] = peak_bed_link
    bigwig_ids = args.bigwig_ids if args.bigwig_ids else [''] * len(pos_bigwigs)
    if len(pos_bigwigs) == len(neg_bigwigs) == len(bigwig_ids):
        for i, (bw_pos, bw_neg, basename) in enumerate(zip(pos_bigwigs, neg_bigwigs, bigwig_ids), start=1):
            basename = basename if basename else right_replace(os.path.basename(bw_pos), '.bw', '')
            bw_pos_link = link_file(bw_pos, os.path.join(outdir, f'{basename}.pos.bw'))
            bw_neg_link = link_file(bw_neg, os.path.join(outdir, f'{basename}.neg.bw'))
            positive_bigwigs[basename] = bw_pos_link
            negative_bigwigs[basename] = bw_neg_link
    return beds, positive_bigwigs, negative_bigwigs, outdir, args.bases


beds, positive_bigwigs, negative_bigwigs, outdir, bases = validate_paths()
names = list(positive_bigwigs.keys())
POS_HANDLERS = {k: pyBigWig.open(v) for k, v in positive_bigwigs.items()}
NEG_HANDLERS = {k: pyBigWig.open(v) for k, v in negative_bigwigs.items()}


class Peak:
    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.peak = f'{chrom}:{start}-{end}:{strand}'


@task(inputs=[], outputs=os.path.join(outdir, 'peaks.union.bed'), kind='create')
def concatenate_bed(inputs, bed):
    cmder.run(f'cat {" ".join([v for v in beds.values()])} > {bed}', msg='Concatenating peak bed files ...')
    return bed


@task(inputs=concatenate_bed, outputs=lambda i: right_replace(i, '.union.bed', '.merged.bed'))
def merge_bed(bed, out):
    tmp = right_replace(bed, '.union.bed', '.tmp.bed')
    cmder.run(f'sort -k1,1 -k2,2n {bed} > {tmp}', msg=f'Sorting {bed} ...')
    # Not output strand ?
    cmder.run(f'bedtools merge -i {tmp} -s -c 6 -o distinct > {out}', msg=f'Merging peaks in {tmp} ...')
    os.unlink(tmp)
    return out


def get_density(chrom, start, end, handler):
    density = handler.values(chrom, start, end + 1)
    return density


@task(inputs=merge_bed, outputs=os.path.join(outdir, 'peaks.density.matrix.tsv'))
def create_matrix(bed, tsv):
    bed = right_replace(tsv, '.density.matrix.tsv', '.merged.bed')
    df = pd.read_csv(bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'strand'])
    dd, headers = [], ['Peak', 'Base']

    for name in names:
        headers.extend([f'{name}:{i}' for i in list(range(-args.bases, 1)) + list(range(1, args.bases + 1))])
    for i, row in enumerate(df.itertuples()):
        densities = []
        peak = Peak(row.chrom, row.start, row.end, row.strand)
        dx = {k: get_density(row.chrom, row.start, row.end, POS_HANDLERS[k] if row.strand == '+' else NEG_HANDLERS[k])
              for k in names}
        dx = pd.DataFrame(dx)
        center = dx.abs().max(axis=1).idxmax() + row.start
        left, right = center - bases, center + bases
        densities.extend([peak.peak, f'{left}-{right}'])
        dx = {k: get_density(row.chrom, left - 1, right - 1, POS_HANDLERS[k] if row.strand == '+' else NEG_HANDLERS[k])
              for k in names}
        maximum = pd.DataFrame(dx).abs().max().max()
        for k in sorted(names):
            densities.extend([abs(v) / maximum for v in dx[k]])
            # densities.extend([abs(v) for v in dx[k]])
        if args.debug and i == 100:
            break
        dd.append(densities)
    dd = pd.DataFrame(dd, columns=headers)
    dd = dd.fillna(0)
    dd = dd[dd.sum(axis=1) > 0]
    dd.to_csv(tsv, index=False, sep='\t', float_format='%.4f')


@task(inputs=create_matrix, outputs=lambda i: i.replace('.density.matrix.tsv', '.heatmap.png'))
def plot_peak(tsv, png):
    df = pd.read_csv(tsv, sep='\t')
    # print(df)
    df = df.drop(columns=['Peak', 'Base'])
    # df = df.drop(columns=[0, 1])
    g = sns.clustermap(df, col_cluster=False, figsize=(12, 8), cmap='Greens', cbar_pos=(0.21, 0.9, 0.78, 0.03),
                       xticklabels=False, yticklabels=False, cbar_kws={"orientation": "horizontal"})
    g.ax_row_dendrogram.set_visible(False)
    g.ax_heatmap.set_xticks(list(range(51, len(names) * 100, 101)))
    g.ax_heatmap.set_xticklabels([name.replace('QKI_', '').replace('.merged.rmDup.r2.sorted.bam.p.sort', '')
                                  for name in names], rotation=90)
    g.savefig(png, dpi=300)
    # g.savefig(right_replace(png, '.png', '.pdf'))
    # g.savefig(right_replace(png, '.png', '.svg'))


def main():
    try:
        flow = Flow('PeakPlot', description=__doc__.strip())
        flow.run(dry=args.dry_run)
        if NEED_TO_DELETE and args.cleanup:
            logger.info('Cleaning up ...')
            for file in need_to_remove:
                cmder.run(f'rm {file}')
            logger.info('Cleaning up complete.')
    finally:
        for k, v in POS_HANDLERS.items():
            v.close()
        for k, v in NEG_HANDLERS.items():
            v.close()


if __name__ == '__main__':
    main()
