# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import subprocess
import shlex
import sys
import os

import click
from . import cli
from ..util import cmd_exists

try:
    from subprocess import DEVNULL  # py3
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')


AWK_TEMPLATE = """\
BEGIN {{
    OFS="\\t";
    print "Enumerating chroms..." > "/dev/stderr";
    i = 1;
    while (getline < "{CHROMSIZES_FILE}") {{
        chrID[$1] = i;
        print $1, i > "/dev/stderr";
        i = i + 1;
    }}
    close("{CHROMSIZES_FILE}");
}}
{{
    if ( !(chrID[${C1}]) || !(chrID[${C2}]) )
        next;
    else if ( (chrID[${C1}] > chrID[${C2}]) || ((chrID[${C1}]==chrID[${C2}]) && (${P1} > ${P2})) )
        print ${C2},${P2},${S2},${C1},${P1},${S1};
    else
        print ${C1},${P1},${S1},${C2},${P2},${S2};
}}"""


@cli.command()
@click.argument(
    "chromsizes_path",
    type=click.Path(exists=True),
    metavar="CHROMSIZES_PATH")
@click.argument(
    "pairs_path",
    type=click.Path(exists=True),
    metavar="PAIRS_PATH")
@click.option(
    "--chrom1", "-c1",
    help="chrom1 field number in the input file (starting from 1)",
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
    "--nproc", "-p",
    help="number of processors",
    type=int,
    default=8)
@click.option(
    "--sort-options",
    help="quoted list of options to `sort`",
    type=str)
@click.option(
    "--out", "-o",
    help="Output gzip file")
def csort(chromsizes_path, pairs_path, chrom1, pos1, strand1, chrom2, pos2,
          strand2, nproc, sort_options, out):
    """
    Sort and index a contact list.
    Arrange the reads of each pair so that all contacts are upper triangular
    with respect to the chromosome ordering given by the chromsizes file, and
    sort contacts by position.

    Requires Unix tools: pigz, sort, bgzip, tabix

    CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
    order.

    PAIRS_PATH : Contacts (i.e. read pairs) text file.

    The output file will have the following properties:

    \b
    - Upper triangular: the read pairs on each row are assigned to side 1 or 2
      in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2),
      according to the desired chromosome order as given by the chromsizes
      file.
    - Rows are lexically sorted by chrom1, pos1, chrom2, pos2. Here, the way
      chromosomes are sorted does not need to respect the desired order.
    - Compressed with bgzip [*]
    - Indexed using Tabix [*] on chrom1 and pos1: `tabix -0 -s1 -b2 -e2`

    [*] Tabix manpage: <http://www.htslib.org/doc/tabix.html>.

    """
    # Check for required Unix tools
    for tool in ['awk', 'sort', 'tabix', 'bgzip']:
        if not cmd_exists(tool):
            print('Command {} not found'.format(tool), file=sys.stderr)
            sys.exit(1)

    # If output path is not given, produce output path by stripping any .txt,
    # .gz or .txt.gz extension from the input path and appending .sorted[.txt].gz
    infile = pairs_path
    if out is None:
        prefix = infile
        ext = '.gz'
        if prefix.endswith('.gz'):
            prefix = op.splitext(prefix)[0]
        if prefix.endswith('.txt'):
            prefix = op.splitext(prefix)[0]
            ext = '.txt.gz'
        outfile = prefix + '.sorted' + ext
    else:
        outfile = out

    # If input file appears gzipped based its extension, read using one of pigz 
    # or gzip for decompression. Otherwise assume uncompressed and use cat.
    ingzip = infile.endswith('.gz')
    if ingzip:
        if cmd_exists('pigz'):
            read_cmd = ['pigz', '-p', str(nproc//2), '-dc',  infile]
        elif cmd_exists('gzip'):
            read_cmd = ['gzip', '-dc', infile]
        else:
            print('No gzip decompressor found.', file=sys.stderr)
            sys.exit(1)
    else:
        read_cmd = ['cat', infile]

    # Generate awk program to reorder paired ends.
    params = {
        'CHROMSIZES_FILE': chromsizes_path,
        'C1': chrom1,
        'P1': pos1,
        'S1': strand1,
        'C2': chrom2,
        'P2': pos2,
        'S2': strand2,
    }
    triu_reorder_cmd = ['awk', AWK_TEMPLATE.format(**params)]

    # Generate sort command
    os.environ['LC_ALL'] = 'C'
    if sort_options is not None:
        sort_options = shlex.split(sort_options)
    elif subprocess.call(['sort', '--parallel=1'],
                         stdin=DEVNULL, 
                         stdout=DEVNULL,
                         stderr=DEVNULL) == 0:
        sort_options = [
            '--parallel={}'.format(nproc//2), 
            '--buffer-size=1G'
        ]
    else:
        sort_options = []
    sort_cmd = ['sort', '-k1,1', '-k2,2n', '-k4,4', '-k5,5n'] + sort_options

    # Run pipeline: re-order reads, sort, then bgzip
    print("Reordering paired end fields and sorting...", file=sys.stderr)
    print("Input: '{}'".format(infile), file=sys.stderr)
    print("Output: '{}'".format(outfile), file=sys.stderr)
    assert infile != outfile
    with open(outfile, 'wb') as f:
        p1 = subprocess.Popen(
            read_cmd,
            stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            triu_reorder_cmd,
            stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(
            sort_cmd,
            stdin=p2.stdout, stdout=subprocess.PIPE)
        p4 = subprocess.Popen(
            ['bgzip', '-c'],
            stdin=p3.stdout, stdout=f)
        p4.communicate()

    # Create tabix index file
    print("Indexing...", file=sys.stderr)
    p5 = subprocess.Popen(['tabix', '-0', '-s1', '-b2', '-e2', outfile])
    p5.communicate()
