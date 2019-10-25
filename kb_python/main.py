import argparse
import logging
import os
import shutil
import sys

from . import __version__
from .config import REFERENCES_MAPPING, TEMP_DIR
from .count import count, count_lamanno
from .ref import download_reference, ref, ref_lamanno
from .utils import get_bustools_version, get_kallisto_version


def display_info():
    kallisto_version = '.'.join(str(i) for i in get_kallisto_version())
    bustools_version = '.'.join(str(i) for i in get_bustools_version())
    info = '''kb_python {}
    kallisto: {}
    bustools: {}
    '''.format(__version__, kallisto_version, bustools_version)
    print(info)


def parse_ref(args):
    if args.d is not None:
        download_reference(args.d, args.i, args.g, overwrite=args.overwrite)
    elif args.lamanno:
        ref_lamanno(
            args.fasta,
            args.gtf,
            args.c,
            args.n,
            args.i,
            args.g,
            args.a,
            args.r,
            overwrite=args.overwrite
        )
    else:
        ref(
            args.fasta,
            args.gtf,
            args.c,
            args.i,
            args.g,
            overwrite=args.overwrite
        )


def parse_count(args):
    if args.lamanno:
        count_lamanno(
            args.i,
            args.g,
            args.c,
            args.n,
            args.x,
            args.o,
            args.fastqs,
            args.w,
            threads=args.t,
            memory=args.m,
            overwrite=args.overwrite,
            loom=args.loom,
            h5ad=args.h5ad,
        )
    else:
        count(
            args.i,
            args.g,
            args.x,
            args.o,
            args.fastqs,
            args.w,
            threads=args.t,
            memory=args.m,
            overwrite=args.overwrite,
            loom=args.loom,
            h5ad=args.h5ad,
        )


COMMAND_TO_FUNCTION = {
    'ref': parse_ref,
    'count': parse_count,
}


def setup_info_args(parser, parent):
    parser_info = parser.add_parser(
        'info',
        description='Display package and citation information',
        help='Display package and citation information',
        parents=[parent],
        add_help=False,
    )
    return parser_info


def setup_ref_args(parser, parent):
    parser_ref = parser.add_parser(
        'ref',
        description='Build a kallisto index and transcript-to-gene mapping',
        help='Build a kallisto index and transcript-to-gene mapping',
        parents=[parent],
    )
    required_ref = parser_ref.add_argument_group('required arguments')
    required_ref.add_argument(
        '-i',
        help='Path to the kallisto index to be constructed',
        type=str,
        required=True
    )
    required_ref.add_argument(
        '-g',
        help='Path to transcript-to-gene mapping to be generated',
        type=str,
        required=True
    )
    required_ref.add_argument(
        '-c',
        help=(
            'Path to the cDNA FASTA to be generated '
            '(not needed when downloading)'
        ),
        type=str,
        required='-d' not in sys.argv
    )
    required_lamanno = parser_ref.add_argument_group(
        'required arguments for --lamanno'
    )
    required_lamanno.add_argument(
        '-n',
        help='Path to the intron FASTA to be generated',
        type=str,
        required='--lamanno' in sys.argv
    )
    required_lamanno.add_argument(
        '-a',
        help='Path to generate cDNA transcripts to be captured',
        type=str,
        required='--lamanno' in sys.argv
    )
    required_lamanno.add_argument(
        '-r',
        help='Path to generate intron transcripts to be captured',
        type=str,
        required='--lamanno' in sys.argv
    )

    parser_ref.add_argument(
        '-d',
        help=(
            'Download a pre-built kallisto index (along with all necessary files) '
            'instead of building it locally.'
        ),
        type=str,
        choices=list(REFERENCES_MAPPING.keys()),
        required=False
    )
    parser_ref.add_argument(
        '--lamanno', help='Prepare files for RNA lamanno', action='store_true'
    )
    parser_ref.add_argument(
        '--overwrite',
        help='Overwrite existing kallisto index',
        action='store_true'
    )
    parser_ref.add_argument(
        'fasta',
        help='Genomic FASTA file',
        type=str,
        nargs=1 if '-d' not in sys.argv else '?'
    )
    parser_ref.add_argument(
        'gtf',
        help='Reference GTF file',
        type=str,
        nargs=1 if '-d' not in sys.argv else '?'
    )
    return parser_ref


def setup_count_args(parser, parent):
    # count
    parser_count = parser.add_parser(
        'count',
        description='Generate count matrices from a set of single-cell FASTQ files',  # noqa
        help='Generate count matrices from a set of single-cell FASTQ files',
        parents=[parent],
    )
    required_count = parser_count.add_argument_group('required arguments')
    required_count.add_argument(
        '-i', help='Path to kallisto index', type=str, required=True
    )
    required_count.add_argument(
        '-g',
        help='Path to transcript-to-gene mapping',
        type=str,
        required=True
    )
    required_count.add_argument(
        '-x', help='Single-cell technology used', type=str, required=True
    )
    required_count.add_argument(
        '-o', help='Path to output directory', type=str, required=True
    )
    parser_count.add_argument(
        '-w',
        help=(
            'Path to file of whitelisted barcodes to correct to. '
            'If not provided and bustools supports the technology, '
            'a pre-packaged whitelist is used. If not, the bustools '
            'whitelist command is used.'
        ),
        type=str
    )
    parser_count.add_argument(
        '-t', help='Number of threads to use (default: 8)', type=int, default=8
    )
    parser_count.add_argument(
        '-m', help='Maximum memory used (default: 4G)', type=str, default='4G'
    )
    required_lamanno = parser_count.add_argument_group(
        'required arguments for --lamanno'
    )
    required_lamanno.add_argument(
        '-c',
        help='Path to cDNA transcripts to be captured',
        type=str,
        required='--lamanno' in sys.argv
    )
    required_lamanno.add_argument(
        '-n',
        help='Path to intron transcripts to be captured',
        type=str,
        required='--lamanno' in sys.argv
    )

    parser_count.add_argument(
        '--lamanno', help='Calculate RNA lamanno', action='store_true'
    )
    parser_count.add_argument(
        '--overwrite',
        help='Overwrite existing output.bus file',
        action='store_true'
    )
    conversion_group = parser_count.add_mutually_exclusive_group()
    conversion_group.add_argument(
        '--loom',
        help='Generate loom file from count matrix',
        action='store_true'
    )
    conversion_group.add_argument(
        '--h5ad',
        help='Generate h5ad file from count matrix',
        action='store_true'
    )
    parser_count.add_argument('fastqs', help='FASTQ files', nargs='+')
    return parser_count


def main():
    # Main parser
    parser = argparse.ArgumentParser(
        description='kb_python {}'.format(__version__)
    )
    subparsers = parser.add_subparsers(
        dest='command',
        title='Where <CMD> can be one of',
        required=True,
    )

    # Add common options to this parent parser
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument(
        '--keep-tmp',
        help='Do not delete the tmp directory',
        action='store_true'
    )
    parent.add_argument(
        '--verbose', help='Print debugging information', action='store_true'
    )

    # Command parsers
    setup_info_args(subparsers, argparse.ArgumentParser(add_help=False))
    parser_ref = setup_ref_args(subparsers, parent)
    parser_count = setup_count_args(subparsers, parent)

    command_to_parser = {
        'ref': parser_ref,
        'count': parser_count,
    }

    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 2:
        if sys.argv[1] == 'info':
            display_info()
        elif sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    logging.basicConfig(
        format='[%(asctime)s] %(levelname)7s %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO,
    )
    logger = logging.getLogger(__name__)
    logger.debug('Printing verbose output')
    logger.debug('Creating tmp directory')
    os.makedirs(TEMP_DIR, exist_ok=True)
    try:
        COMMAND_TO_FUNCTION[args.command](args)
    finally:
        # Always clean temp dir
        if not args.keep_tmp:
            shutil.rmtree(TEMP_DIR, ignore_errors=True)
