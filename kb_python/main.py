import argparse
import sys

from . import __version__
from .count import count
from .ref import ref


def parse_ref(args):
    ref(args.fasta, args.gtf, args.i, args.g, overwrite=args.overwrite)


def parse_count(args):
    count(
        args.i,
        args.g,
        args.x,
        args.o,
        args.fastqs,
        args.w,
        threads=args.t,
        memory=args.m,
        keep_temp=args.keep_tmp,
        overwrite=args.overwrite,
        loom=args.loom,
        h5ad=args.h5ad,
    )


COMMAND_TO_FUNCTION = {
    'ref': parse_ref,
    'count': parse_count,
}


def setup_ref_args(parser):
    parser_ref = parser.add_parser(
        'ref',
        description='Build a kallisto index and transcript-to-gene mapping',
        help='Build a kallisto index and transcript-to-gene mapping',
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
    parser_ref.add_argument(
        '--overwrite',
        help='Overwrite existing kallisto index',
        action='store_true'
    )
    parser_ref.add_argument('fasta', help='Reference FASTA file', type=str)
    parser_ref.add_argument('gtf', help='Reference GTF file', type=str)
    return parser_ref


def setup_count_args(parser):
    # count
    parser_count = parser.add_parser(
        'count',
        description='Generate count matrices from a set of single-cell FASTQ files',  # noqa
        help='Generate count matrices from a set of single-cell FASTQ files',
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
    parser_count.add_argument(
        '--keep-tmp',
        help='Do not delete the tmp directory',
        action='store_true'
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
    parser = argparse.ArgumentParser(
        description='kb_python {}'.format(__version__),
    )
    subparsers = parser.add_subparsers(
        dest='command',
        title='Where <CMD> can be one of',
        required=True,
    )
    parser_ref = setup_ref_args(subparsers)
    parser_count = setup_count_args(subparsers)

    command_to_parser = {
        'ref': parser_ref,
        'count': parser_count,
    }

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    COMMAND_TO_FUNCTION[args.command](args)
