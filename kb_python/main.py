import argparse
import logging
import os
import shutil
import sys
import textwrap

from . import __version__
from .config import PACKAGE_PATH, REFERENCES_MAPPING, TECHNOLOGIES, TEMP_DIR
from .constants import INFO_FILENAME
from .count import count, count_velocity
from .ref import download_reference, ref, ref_lamanno
from .utils import get_bustools_version, get_kallisto_version


def display_info():
    """Displays kb, kallisto and bustools version + citation information, along
    with a brief description and examples.
    """
    kallisto_version = '.'.join(str(i) for i in get_kallisto_version())
    bustools_version = '.'.join(str(i) for i in get_bustools_version())
    info = '''kb_python {}
    kallisto: {}
    bustools: {}
    '''.format(__version__, kallisto_version, bustools_version)
    with open(os.path.join(PACKAGE_PATH, INFO_FILENAME), 'r') as f:
        print(
            '{}\n{}'.format(
                info, '\n'.join([
                    line.strip()
                    if line.startswith('(') else textwrap.fill(line, width=80)
                    for line in f.readlines()
                ])
            )
        )
    sys.exit(1)


def display_technologies():
    """Displays a list of supported technologies along with whether kb provides
    a whitelist for that technology and the FASTQ argument order for kb count.
    """
    headers = [
        'name', 'whitelist provided', 'barcode (file #, start, stop)',
        'umi (file #, start, stop)', 'read file #'
    ]
    rows = [headers]

    print('List of supported single-cell technologies\n')
    for t in TECHNOLOGIES:
        row = [
            t.name,
            'yes' if t.whitelist_archive else '',
            ' '.join(str(tup) for tup in t.barcode_positions),
            ' '.join(str(tup) for tup in t.umi_positions),
            str(t.reads_file),
        ]
        rows.append(row)

    max_lens = []
    for i in range(len(headers)):
        max_lens.append(len(headers[i]))
        for row in rows[1:]:
            max_lens[i] = max(max_lens[i], len(row[i]))

    rows.insert(1, ['-' * l for l in max_lens])
    for row in rows:
        for col, l in zip(row, max_lens):
            print(col.ljust(l + 4), end='')
        print()
    sys.exit(1)


def parse_ref(args):
    """Parser for the `ref` command.

    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    if args.d is not None:
        # Options that are files.
        options = ['i', 'g', 'c1', 'c2']
        files = {
            option: getattr(args, option)
            for option in options
            if getattr(args, option) is not None
        }
        reference = REFERENCES_MAPPING[args.d]
        download_reference(reference, files, overwrite=args.overwrite)
    elif args.lamanno:
        ref_lamanno(
            args.fasta,
            args.gtf,
            args.f1,
            args.f2,
            args.i,
            args.g,
            args.c1,
            args.c2,
            overwrite=args.overwrite
        )
    else:
        ref(
            args.fasta,
            args.gtf,
            args.f1,
            args.i,
            args.g,
            overwrite=args.overwrite
        )


def parse_count(args):
    """Parser for the `count` command.

    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    if args.lamanno or args.nucleus:
        count_velocity(
            args.i,
            args.g,
            args.c1,
            args.c2,
            args.x,
            args.o,
            args.fastqs,
            args.w,
            filter=args.filter,
            threads=args.t,
            memory=args.m,
            overwrite=args.overwrite,
            loom=args.loom,
            h5ad=args.h5ad,
            nucleus=args.nucleus,
        )
    else:
        count(
            args.i,
            args.g,
            args.x,
            args.o,
            args.fastqs,
            args.w,
            filter=args.filter,
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
    """Helper function to set up a subparser for the `info` command.

    :param parser: argparse parser to add the `info` command to
    :type args: argparse.ArgumentParser
    :param parent: argparse parser parent of the newly added subcommand.
                   used to inherit shared commands/flags
    :type args: argparse.ArgumentParser
    :return: the newly added parser
    :rtype: argparse.ArgumentParser
    """
    parser_info = parser.add_parser(
        'info',
        description='Display package and citation information',
        help='Display package and citation information',
        parents=[parent],
        add_help=False,
    )
    return parser_info


def setup_ref_args(parser, parent):
    """Helper function to set up a subparser for the `ref` command.

    :param parser: argparse parser to add the `ref` command to
    :type args: argparse.ArgumentParser
    :param parent: argparse parser parent of the newly added subcommand.
                   used to inherit shared commands/flags
    :type args: argparse.ArgumentParser
    :return: the newly added parser
    :rtype: argparse.ArgumentParser
    """
    parser_ref = parser.add_parser(
        'ref',
        description='Build a kallisto index and transcript-to-gene mapping',
        help='Build a kallisto index and transcript-to-gene mapping',
        parents=[parent],
    )
    parser_ref._actions[0].help = parser_ref._actions[0].help.capitalize()

    required_ref = parser_ref.add_argument_group('required arguments')
    required_ref.add_argument(
        '-i',
        metavar='INDEX',
        help='Path to the kallisto index to be constructed',
        type=str,
        required=True
    )
    required_ref.add_argument(
        '-g',
        metavar='T2G',
        help='Path to transcript-to-gene mapping to be generated',
        type=str,
        required=True
    )
    required_ref.add_argument(
        '-f1',
        metavar='FASTA',
        help=('[Optional with -d] Path to the cDNA FASTA to be generated '),
        type=str,
        required='-d' not in sys.argv
    )
    required_lamanno = parser_ref.add_argument_group(
        'required arguments for --lamanno'
    )
    required_lamanno.add_argument(
        '-f2',
        metavar='FASTA',
        help='Path to the intron FASTA to be generated',
        type=str,
        required='--lamanno' in sys.argv
    )
    required_lamanno.add_argument(
        '-c1',
        metavar='T2C',
        help='Path to generate cDNA transcripts-to-capture',
        type=str,
        required='--lamanno' in sys.argv
    )
    required_lamanno.add_argument(
        '-c2',
        metavar='T2C',
        help='Path to generate intron transcripts-to-capture',
        type=str,
        required='--lamanno' in sys.argv
    )

    parser_ref.add_argument(
        '-d',
        help=(
            'Download a pre-built kallisto index (along with all necessary files) '
            'instead of building it locally'
        ),
        type=str,
        choices=list(REFERENCES_MAPPING.keys()),
        required=False
    )
    parser_ref.add_argument(
        '--lamanno',
        help=(
            'Prepare files for RNA velocity based on '
            'La Manno et al. 2018 logic'
        ),
        action='store_true'
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
        nargs=None if '-d' not in sys.argv else '?'
    )
    parser_ref.add_argument(
        'gtf',
        help='Reference GTF file',
        type=str,
        nargs=None if '-d' not in sys.argv else '?'
    )
    return parser_ref


def setup_count_args(parser, parent):
    """Helper function to set up a subparser for the `count` command.

    :param parser: argparse parser to add the `count` command to
    :type args: argparse.ArgumentParser
    :param parent: argparse parser parent of the newly added subcommand.
                   used to inherit shared commands/flags
    :type args: argparse.ArgumentParser
    :return: the newly added parser
    :rtype: argparse.ArgumentParser
    """
    # count
    parser_count = parser.add_parser(
        'count',
        description=('Generate count matrices from a set of single-cell FASTQ files. '
                     'Run `kb --list` to view single-cell technology information.'),  # noqa
        help='Generate count matrices from a set of single-cell FASTQ files',
        parents=[parent],
    )
    parser_count._actions[0].help = parser_count._actions[0].help.capitalize()

    required_count = parser_count.add_argument_group('required arguments')
    required_count.add_argument(
        '-i',
        metavar='INDEX',
        help='Path to kallisto index',
        type=str,
        required=True
    )
    required_count.add_argument(
        '-g',
        metavar='T2G',
        help='Path to transcript-to-gene mapping',
        type=str,
        required=True
    )
    required_count.add_argument(
        '-x',
        metavar='TECHNOLOGY',
        help='Single-cell technology used (`kb --list` to view)',
        type=str,
        required=True
    )
    parser_count.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    parser_count.add_argument(
        '-w',
        metavar='WHITELIST',
        help=(
            'Path to file of whitelisted barcodes to correct to. '
            'If not provided and bustools supports the technology, '
            'a pre-packaged whitelist is used. If not, the bustools '
            'whitelist command is used. (`kb --list` to view whitelists)'
        ),
        type=str
    )
    parser_count.add_argument(
        '-t',
        metavar='THREADS',
        help='Number of threads to use (default: 8)',
        type=int,
        default=8
    )
    parser_count.add_argument(
        '-m',
        metavar='MEMORY',
        help='Maximum memory used (default: 4G)',
        type=str,
        default='4G'
    )
    required_lamanno = parser_count.add_argument_group(
        'required arguments for --lamanno and --nucleus'
    )
    required_lamanno.add_argument(
        '-c1',
        metavar='T2C',
        help='Path to cDNA transcripts-to-capture',
        type=str,
        required=any(arg in sys.argv for arg in ('--lamanno', '--nucleus'))
    )
    required_lamanno.add_argument(
        '-c2',
        metavar='T2C',
        help='Path to intron transcripts-to-captured',
        type=str,
        required=any(arg in sys.argv for arg in ('--lamanno', '--nucleus'))
    )
    parser_count.add_argument(
        '--overwrite',
        help='Overwrite existing output.bus file',
        action='store_true'
    )

    velocity_group = parser_count.add_mutually_exclusive_group()
    velocity_group.add_argument(
        '--lamanno',
        help='Calculate RNA velocity based on La Manno et al. 2018 logic',
        action='store_true'
    )
    velocity_group.add_argument(
        '--nucleus',
        help='Calculate RNA velocity on single-nucleus RNA-seq reads',
        action='store_true'
    )
    parser_count.add_argument(
        '--filter',
        help='Produce a filtered gene count matrix (default: bustools)',
        type=str,
        const='bustools',
        nargs='?',
        choices=['bustools']
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
    """Command-line entrypoint.
    """
    # Main parser
    parser = argparse.ArgumentParser(
        description='kb_python {}'.format(__version__)
    )
    parser._actions[0].help = parser._actions[0].help.capitalize()
    parser.add_argument(
        '--list',
        help='Display list of supported single-cell technologies',
        action='store_true'
    )
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='<CMD>',
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
    if 'info' in sys.argv:
        display_info()
    elif '--list' in sys.argv:
        display_technologies()

    # Show help when no arguments are given
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

    logging.basicConfig(
        format='[%(asctime)s] %(levelname)7s %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO,
    )
    logger = logging.getLogger(__name__)
    logger.debug('Printing verbose output')
    logger.debug('Creating tmp directory')
    os.makedirs(TEMP_DIR, exist_ok=True)
    try:
        logger.debug(args)
        COMMAND_TO_FUNCTION[args.command](args)
    except Exception:
        logger.exception('An exception occurred')
    finally:
        # Always clean temp dir
        if not args.keep_tmp:
            logger.debug('Removing tmp directory')
            shutil.rmtree(TEMP_DIR, ignore_errors=True)
