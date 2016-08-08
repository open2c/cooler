#!/usr/bin/env python
from __future__ import division, print_function
import os.path as op
import subprocess
import argparse
import sys
import os


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
    if ( (chrID[${C1}] > chrID[${C2}]) || ((chrID[${C1}]==chrID[${C2}]) && (${P1} > ${P2})) )
        print ${C2},${P2},${S2},${C1},${P1},${S1};
    else
        print ${C1},${P1},${S1},${C2},${P2},${S2};
}}"""


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Sort contacts by position and order the reads of each pair so that all "
                    "contacts are upper triangular with respect to the chromosome ordering "
                    "given by the chromsizes file.")
    parser.add_argument(
        "chromsizes",
        help="UCSC-like chromsizes file, with chromosomes in desired order",
        metavar="CHROMSIZES_PATH")
    parser.add_argument(
        "pairs",
        help="Contact list file",
        metavar="PAIRS_PATH")
    parser.add_argument(
        "--chrom1", "-c1",
        help="chrom1 field number",
        default=1)
    parser.add_argument(
        "--pos1", "-p1",
        help="pos1 field number",
        default=2)
    parser.add_argument(
        "--strand1", "-s1",
        help="strand1 field number",
        default=3)
    parser.add_argument(
        "--chrom2", "-c2",
        help="chrom2 field number",
        default=4)
    parser.add_argument(
        "--pos2", "-p2",
        help="pos2 field number",
        default=5)
    parser.add_argument(
        "--strand2", "-s2",
        help="strand2 field number",
        default=6)
    parser.add_argument(
        "--nproc", "-p",
        help="number of processors",
        type=int,
        default=8)
    parser.add_argument(
        "--out", "-o",
        help="Output gzip file")

    args = vars(parser.parse_args())

    chromsizes = op.realpath(args['chromsizes'])
    infile = op.realpath(args['pairs'])
    if args['out'] is None:
        outfile = infile.replace('.txt.gz', 'sorted.txt.gz')
    else:
        outfile = args['out']

    params = {
        'CHROMSIZES_FILE': chromsizes,
        'C1': args['chrom1'],
        'P1': args['pos1'],
        'S1': args['strand1'],
        'C2': args['chrom2'],
        'P2': args['pos2'],
        'S2': args['strand2'],
    }
    triu_reorder = AWK_TEMPLATE.format(**params)

    nproc = args['nproc']
    os.environ['LC_ALL'] = 'C'

    # Re-order reads, sort, then bgzip
    with open(outfile, 'wb') as f:
        p1 = subprocess.Popen(
            ['pigz', '-p', nproc//2, '-dc',  infile],
            stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ['awk', triu_reorder],
            stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(
            ['sort', '-S1G', '--parallel={}'.format(nproc//2), '-k1,1', '-k2,2n', '-k4,4', '-k5,5n'],
            stdin=p2.stdout, stdout=subprocess.PIPE)
        p4 = subprocess.Popen(
            ['bgzip', '-c'],
            stdin=p3.stdout, stdout=f)
        p4.communicate()

    # Create tabix index file
    p5 = subprocess.Popen(['tabix', '-0', '-s1', '-b2', '-e2', outfile])
    p5.communicate()
