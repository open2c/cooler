# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import subprocess
import sys
import os

import click
from . import cli


AWK_TEMPLATE = """\
BEGIN {{
    OFS="\\t";
    i = 0;
    while (getline < "{CHROMSIZES_FILE}") {{
        chrID[${C1}] = i;
        i = i + 1;
    }}
    close("{CHROMSIZES_FILE}");
}}
{{
    if ( (chrID[${C1}] < chrID[${C2}]) || ((chrID[${C1}]==chrID[${C2}]) && (${P1} > ${P2})) )
        print ${C2},${P2},${S2},${C1},${P1},${S1};
    else
        print ${C1},${P1},${S1},${C2},${P2},${S2};
}}"""


@cli.command()
@click.argument(
    "chromsizes_path",
    metavar="CHROMSIZES_PATH")
@click.argument(
    "pairs_path",
    metavar="PAIRS_PATH")
@click.option(
    "--chrom1", "-c1",
    help="chrom1 field number",
    default=1)
@click.option(
    "--pos1", "-p1",
    help="pos1 field number",
    default=2)
@click.option(
    "--strand1", "-s1",
    help="strand1 field number",
    default=3)
@click.option(
    "--chrom2", "-c2",
    help="chrom2 field number",
    default=4)
@click.option(
    "--pos2", "-p2",
    help="pos2 field number",
    default=5)
@click.option(
    "--strand2", "-s2",
    help="strand2 field number",
    default=6)
@click.option(
    "--ncpu", "-p",
    help="number of processors",
    type=int,
    default=8)
@click.option(
    "--out", "-o",
    help="Output gzip file")
def csort(chromsizes_path, pairs_path, chrom1, pos1, strand1, chrom2, pos2, strand2, out):
    """
    Sort contacts by position and order the reads of each pair so that all 
    contacts are upper triangular with respect to the chromosome ordering 
    given by the chromsizes file.

    CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired order.

    PAIRS_PATH : Read pairs (i.e. contacts) text file.

    """
    chromsizes_path = op.realpath(chromsizes_path)
    infile = op.realpath(pairs_path)
    if out is None:
        outfile = infile.replace('.txt.gz', 'sorted.txt.gz')
    else:
        outfile = out

    params = {
        'CHROMSIZES_FILE': chromsizes_path,
        'C1': chrom1,
        'P1': pos1,
        'S1': strand1,
        'C2': chrom2,
        'P2': pos2,
        'S2': strand2,
    }
    triu_reorder = AWK_TEMPLATE.format(**params)

    os.environ['LC_ALL'] = 'C'

    # Re-order reads, sort, then bgzip
    with open(outfile, 'wb') as f:
        p1 = subprocess.Popen(
            ['pigz', '-p', ncpu//2, '-dc',  infile],
            stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ['awk', triu_reorder],
            stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(
            ['sort', '-S1G', '--parallel={}'.format(ncpu//2), '-k1,1', '-k2,2n', '-k4,4', '-k5,5n'],
            stdin=p2.stdout, stdout=subprocess.PIPE)
        p4 = subprocess.Popen(
            ['bgzip', '-c'],
            stdin=p3.stdout, stdout=f)
        p4.communicate()

    # Create tabix index file
    p5 = subprocess.Popen(['tabix', '-0', '-s1', '-b2', '-e2', outfile])
    p5.communicate()
