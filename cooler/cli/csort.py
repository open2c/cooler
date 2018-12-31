# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import subprocess
import shlex
import sys
import os

from . import cli, get_logger
import click


from ..util import cmd_exists

try:
    from subprocess import DEVNULL  # py3
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')


logger = get_logger(__name__)


SORT_POS = 'sort -k{C1},{C1} -k{P1},{P1}n -k{C2},{C2} -k{P2},{P2}n'

SORT_BLK = 'sort -k{C1},{C1} -k{C2},{C2} -k{P1},{P1}n -k{P2},{P2}n'

INDEX_TBX = 'tabix -f -s{C1} -b{P1} -e{P1}'

INDEX_PX2 = 'pairix -f -s{C1} -d{C2} -b{P1} -e{P1} -u{P2} -v{P2}'

FLIP_TEMPLATE = \
"""import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

instream = getattr(sys.stdin, "buffer", sys.stdin)
outstream = getattr(sys.stdout, "buffer", sys.stdout)

with open("{chromosomes_path}", "rb") as f:
    chrIDs = {{}}
    for i, line in enumerate(f, 1):
        chrom = line.split(b"\\t")[0].strip()
        if chrom:
            chrIDs[chrom] = i

for line in instream:
    if line.startswith(b"{comment_char}"):
        continue
    parts = line.strip().split(b"{sep}")
    chrom1, chrom2 = parts[{c1}], parts[{c2}]
    if chrom1 in chrIDs and chrom2 in chrIDs:
        cid1, cid2 = chrIDs[chrom1], chrIDs[chrom2]
        pos1, pos2 = int(parts[{p1}]), int(parts[{p2}])
        if (cid1 > cid2) or ((cid1 == cid2) and (pos1 > pos2)):
            parts[{c1}], parts[{c2}] = parts[{c2}], parts[{c1}]
            parts[{p1}], parts[{p2}] = parts[{p2}], parts[{p1}]
        outstream.write(b"\\t".join(parts) + b"\\n")
"""


def _has_parallel_sort():
    return subprocess.call(['sort', '--parallel=1'],
                           stdin=DEVNULL,
                           stdout=DEVNULL,
                           stderr=DEVNULL) == 0


def _validate_fieldnum(ctx, param, value):
    if value <= 0:
        raise click.BadParameter('Field numbers are one-based')
    return value


def make_read_command(infile):
    """
    If input file appears gzipped based its extension, read using one of pigz
    or gzip for decompression. Otherwise assume uncompressed and use cat.

    Note
    ----
    Gzip decompression can't actually be parallelized, but pigz will create
    three extra threads that may help.
    See <https://github.com/madler/pigz/issues/36>.

    """
    ingzip = infile.endswith('.gz')
    if ingzip:
        if cmd_exists('pigz'):
            read_cmd = ['pigz', '-dc', infile]
        elif cmd_exists('gzip'):
            read_cmd = ['gzip', '-dc', infile]
        else:
            print('No gzip decompressor found.', file=sys.stderr)
            sys.exit(1)
    else:
        read_cmd = ['cat', infile]

    return read_cmd


def make_flip_command(chromosomes_path, sep, comment_char, fields):
    with open(chromosomes_path, 'rt') as f:
        logger.info("Enumerating requested chromosomes...")
        for i, line in enumerate(f, 1):
            chrom = line.split('\t')[0].strip()
            if chrom:
                logger.info(chrom + '\t' + str(i))

    # zero-based column indices
    c1 = fields['C1'] - 1
    c2 = fields['C2'] - 1
    p1 = fields['P1'] - 1
    p2 = fields['P2'] - 1
    flip_code = FLIP_TEMPLATE.format(
        chromosomes_path=chromosomes_path,
        sep=sep,
        comment_char=comment_char,
        c1=c1, c2=c2, p1=p1, p2=p2)

    return [sys.executable, '-c', flip_code]


def make_sort_command(index, fields, sort_options):
    os.environ['LC_ALL'] = 'C'
    if index == 'tabix':
        sort_cmd = shlex.split(SORT_POS.format(**fields))
    elif index == 'pairix':
        sort_cmd = shlex.split(SORT_BLK.format(**fields))
    sort_cmd += sort_options
    return sort_cmd


def make_index_command(index, fields, zero_based, outfile):
    if index == 'tabix':
        index_cmd = shlex.split(INDEX_TBX.format(**fields))
    elif index == 'pairix':
        index_cmd = shlex.split(INDEX_PX2.format(**fields))
    if zero_based:
        index_cmd += ['-0']
    return index_cmd + [outfile]


@cli.command()
@click.argument(
    "pairs_path",
    type=click.Path(exists=True, allow_dash=True),
    metavar="PAIRS_PATH")
@click.argument(
    "chromosomes_path",
    type=click.Path(exists=True),
    metavar="CHROMOSOMES_PATH")
@click.option(
    "--chrom1", "-c1",
    required=True,
    help="chrom1 field number in the input file (starting from 1)",
    type=int,
    callback=_validate_fieldnum)
@click.option(
    "--chrom2", "-c2",
    required=True,
    help="chrom2 field number",
    type=int,
    callback=_validate_fieldnum)
@click.option(
    "--pos1", "-p1",
    required=True,
    help="pos1 field number",
    type=int,
    callback=_validate_fieldnum)
@click.option(
    "--pos2", "-p2",
    required=True,
    help="pos2 field number",
    type=int,
    callback=_validate_fieldnum)
@click.option(
    "--index", "-i",
    help="Select the preset sort and indexing options",
    type=click.Choice(['tabix', 'pairix']),
    default='pairix',
    show_default=True)
@click.option(
    "--flip-only",
    help="Only flip mates; no sorting or indexing. Write to stdout.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--nproc", "-p",
    help="Number of processors",
    type=int,
    default=8,
    show_default=True)
@click.option(
    "--zero-based", "-0",
    help="Read positions are zero-based",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--sep",
    help="Data delimiter in the input file",
    type=str,
    default=r"\t",
    show_default=True)
@click.option(
    "--comment-char",
    help="Comment character to skip header",
    type=str,
    default='#',
    show_default=True)
@click.option(
    "--sort-options",
    help="Quoted list of additional options to `sort` command",
    type=str)
@click.option(
    "--out", "-o",
    help="Output gzip file")
@click.option(
    "--strand1", "-s1",
    help="strand1 field number (deprecated)",
    type=int)
@click.option(
    "--strand2", "-s2",
    help="strand2 field number (deprecated)",
    type=int)
def csort(pairs_path, chromosomes_path, index, chrom1, chrom2, pos1, pos2,
          flip_only, nproc, zero_based, sep, comment_char, sort_options, out,
          **kwargs):
    """
    Sort and index a contact list.

    Order the mates of each pair record so that all contacts are upper
    triangular with respect to the chromosome ordering given by the chromosomes
    file, sort contacts by genomic location, and index the resulting file.

    PAIRS_PATH : Contacts (i.e. read pairs) text file, optionally compressed.

    CHROMOSOMES_PATH : File listing desired chromosomes in the desired order.
    May be tab-delimited, e.g. a UCSC-style chromsizes file. Contacts mapping to
    other chromosomes will be discarded.

    **Notes**

    \b
    - csort can also be used to sort and index a text representation of
      a contact *matrix* in bedGraph-like format. In this case, substitute
      `pos1` and `pos2` with `start1` and `start2`, respectively.
    - Requires Unix tools: sort, bgzip + tabix or pairix.

    If indexing with Tabix, the output file will have the following properties:

    \b
    - Upper triangular: the read pairs on each row are assigned to side 1 or 2
      in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2)
    - Rows are lexicographically sorted by chrom1, pos1, chrom2, pos2;
      i.e. "positionally sorted"
    - Compressed with bgzip [*]
    - Indexed using Tabix [*] on chrom1 and pos1.

    If indexing with Pairix, the output file will have the following properties:

    \b
    - Upper triangular: the read pairs on each row are assigned to side 1 or 2
      in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2)
    - Rows are lexicographically sorted by chrom1, chrom2, pos1, pos2; i.e.
      "block sorted"
    - Compressed with bgzip [*]
    - Indexed using Pairix [+] on chrom1, chrom2 and pos1.

    \b
    [*] Tabix manpage: <http://www.htslib.org/doc/tabix.html>.
    [+] Pairix on Github: <https://github.com/4dn-dcic/pairix>

    """
    if os.name == 'nt':
        raise click.Abort(
            '"cooler csort" does not work on Windows. To ingest unsorted pairs '
            'data, see the "cooler cload pairs" command.')

    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    # Check for required Unix tools
    for tool in ['sort', 'bgzip'] + [index]:
        if not cmd_exists(tool):
            print('Command {} not found'.format(tool), file=sys.stderr)
            sys.exit(1)

    # If output path is not given, produce output path by stripping any .txt,
    # .gz or .txt.gz extension from the input path and appending .sorted[.txt].gz
    infile = pairs_path
    if out is None:
        if infile == '-' and not flip_only:
            logger.error("Output name required when input is stdin")
            raise click.Abort
        prefix = infile
        ext = '.gz'
        if prefix.endswith('.gz'):
            prefix = op.splitext(prefix)[0]
        if prefix.endswith('.txt'):
            prefix = op.splitext(prefix)[0]
            ext = '.txt.gz'
        if index == 'pairix':
            sort_style = '.blksrt'
        else:
            sort_style = '.possrt'
        outfile = prefix + sort_style + ext
    else:
        outfile = out

    # Parse extra sort options and determine if sort supports --parallel option
    if sort_options is not None:
        sort_options = shlex.split(sort_options)
    elif _has_parallel_sort():
        sort_options = [
            '--parallel={}'.format(nproc),
            '--buffer-size=50%'
        ]
    else:
        sort_options = []

    # 1-based column numbers
    fields = {
        'C1': chrom1,
        'P1': pos1,
        'C2': chrom2,
        'P2': pos2,
    }

    # build commands
    read_cmd = make_read_command(infile)
    flip_cmd = make_flip_command(chromosomes_path, sep, comment_char, fields)

    if flip_only:
        # run pipeline
        logger.info("Reordering pair mates...")
        pipeline = []

        logger.debug(' '.join(read_cmd))
        pipeline.append(
            subprocess.Popen(
                read_cmd,
                stdin=sys.stdin if infile == '-' else None,
                stdout=subprocess.PIPE)
        )

        logger.debug(' '.join(flip_cmd))
        pipeline.append(
            subprocess.Popen(
                flip_cmd,
                stdin=pipeline[-1].stdout,
                stdout=sys.stdout)
        )
        for p in pipeline[::-1]:
            p.communicate()
            if p.returncode != 0:
                sys.exit(1)
    else:

        sort_cmd = make_sort_command(index, fields, sort_options)
        write_cmd = ['bgzip', '-c']
        index_cmd = make_index_command(index, fields, zero_based, outfile)

        # run pipeline
        logger.info("Input: '{}'".format(infile))
        logger.info("Output: '{}'".format(outfile))
        assert infile != outfile

        with open(outfile, 'wb') as fout:

            pipeline = []

            logger.debug(' '.join(read_cmd))
            pipeline.append(
                subprocess.Popen(
                    read_cmd,
                    stdin=sys.stdin if infile == '-' else None,
                    stdout=subprocess.PIPE)
            )

            logger.info("Reordering pair mates and sorting pair records...")
            logger.debug(' '.join(flip_cmd))
            pipeline.append(
                subprocess.Popen(
                    flip_cmd,
                    stdin=pipeline[-1].stdout,
                    stdout=subprocess.PIPE)
            )

            if index == 'pairix':
                logger.info("Sort order: block (chrom1, chrom2, pos1, pos2)")
            else:
                logger.info("Sort order: positional (chrom1, pos1, chrom2, pos2)")
            logger.info(' '.join(sort_cmd))
            pipeline.append(
                subprocess.Popen(
                    sort_cmd,
                    stdin=pipeline[-1].stdout,
                    stdout=subprocess.PIPE)
            )

            logger.debug(' '.join(write_cmd))
            pipeline.append(
                subprocess.Popen(
                    write_cmd,
                    stdin=pipeline[-1].stdout,
                    stdout=fout)
            )

            for p in pipeline[::-1]:
                p.communicate()

                if p.returncode != 0:
                    logger.error(' '.join(p.args))
                    sys.exit(1)

        # Create index file
        logger.info("Indexing...")
        logger.info("Indexer: {}".format(index))
        logger.info(' '.join(index_cmd))
        p = subprocess.Popen(index_cmd)
        p.communicate()
        if p.returncode != 0:
            sys.exit(1)
