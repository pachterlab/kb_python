import argparse
import glob
import logging
import os
import sys
import textwrap
import warnings

from . import __version__
from .config import (
    get_bustools_binary_path,
    get_kallisto_binary_path,
    is_dry,
    no_validate,
    PACKAGE_PATH,
    REFERENCES_MAPPING,
    set_dry,
    set_bustools_binary_path,
    set_kallisto_binary_path,
    TECHNOLOGIES,
    TEMP_DIR,
)
from .constants import INFO_FILENAME
from .count import count, count_smartseq, count_velocity
from .logging import logger
from .ref import download_reference, ref, ref_kite, ref_lamanno
from .utils import (
    get_bustools_version,
    get_kallisto_version,
    make_directory,
    open_as_text,
    remove_directory,
)


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
    headers = ['name', 'whitelist', 'barcode', 'umi', 'cDNA']
    rows = [headers]

    print('List of supported single-cell technologies\n')
    print('Positions syntax: `input file index, start position, end position`')
    print('When start & end positions are None, refers to the entire file')
    print(
        'Custom technologies may be defined by providing a kallisto-supported '
        'technology string\n(see https://pachterlab.github.io/kallisto/manual)\n'
    )
    for t in TECHNOLOGIES:
        chem = t.chemistry
        row = [
            t.name,
            'yes' if chem.has_whitelist else '',
            ' '.join(str(_def) for _def in chem.cell_barcode_parser)
            if chem.has_cell_barcode else '',
            ' '.join(str(_def)
                     for _def in chem.umi_parser) if chem.has_umi else '',
            ' '.join(str(_def) for _def in chem.cdna_parser),
        ]
        rows.append(row)

    max_lens = []
    for i in range(len(headers)):
        max_lens.append(len(headers[i]))
        for row in rows[1:]:
            max_lens[i] = max(max_lens[i], len(row[i]))

    rows.insert(1, ['-' * l for l in max_lens])  # noqa
    for row in rows:
        for col, l in zip(row, max_lens):
            print(col.ljust(l + 4), end='')
        print()
    sys.exit(1)


def parse_ref(parser, args, temp_dir='tmp'):
    """Parser for the `ref` command.

    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    if args.k is not None:
        if args.k < 0 or not args.k % 2:
            parser.error('K-mer length must be a positive odd integer.')
    if args.fasta:
        args.fasta = args.fasta.split(',')
    if args.gtf:
        args.gtf = args.gtf.split(',')
    if (args.fasta and args.gtf) and len(args.fasta) != len(args.gtf):
        parser.error(
            'There must be the same number of FASTAs as there are GTFs.'
        )

    if args.d is not None:
        # Options that are files.
        options = ['i', 'g', 'c1', 'c2']
        files = {
            option: getattr(args, option)
            for option in options
            if getattr(args, option) is not None
        }
        reference = REFERENCES_MAPPING[args.d]
        download_reference(
            reference, files, overwrite=args.overwrite, temp_dir=temp_dir
        )
    elif args.workflow in {'lamanno', 'nucleus'}:
        ref_lamanno(
            args.fasta,
            args.gtf,
            args.f1,
            args.f2,
            args.i,
            args.g,
            args.c1,
            args.c2,
            n=args.n,
            k=args.k,
            flank=args.flank,
            overwrite=args.overwrite,
            temp_dir=temp_dir
        )
    else:
        # Report extraneous options
        velocity_only = ['f2', 'c1', 'c2', 'flank']
        for arg in velocity_only:
            if getattr(args, arg):
                parser.error(
                    f'Option `{arg}` is not supported for workflow `{args.workflow}`'
                )

        if args.workflow == 'kite':
            ref_kite(
                args.feature,
                args.f1,
                args.i,
                args.g,
                n=args.n,
                k=args.k,
                no_mismatches=args.no_mismatches,
                overwrite=args.overwrite,
                temp_dir=temp_dir
            )
        else:
            ref(
                args.fasta,
                args.gtf,
                args.f1,
                args.i,
                args.g,
                n=args.n,
                k=args.k,
                overwrite=args.overwrite,
                temp_dir=temp_dir
            )


def parse_count(parser, args, temp_dir='tmp'):
    """Parser for the `count` command.

    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    if args.report:
        logger.warning((
            'Using `--report` may cause `kb` to exceed maximum memory specified '
            'and crash for large count matrices.'
        ))

    args.i = args.i.split(',')
    if len(args.i) > 1:
        logger.warning((
            'Multiple indices were provided. Aligning to split indices is currently '
            'EXPERIMENTAL and results in loss of reads. It is recommended to '
            'use a single index until this feature is fully supported. Use at '
            'your own risk!'
        ))

    if args.w and args.w.lower() == 'none':
        args.w = None

    if args.workflow in {'lamanno', 'nucleus'}:
        # Smartseq can not be used with lamanno or nucleus.
        if args.x.upper() == 'SMARTSEQ':
            parser.error(
                f'Technology `{args.x}` can not be used with workflow {args.workflow}.'
            )

        count_velocity(
            args.i,
            args.g,
            args.c1,
            args.c2,
            args.x,
            args.o,
            args.fastqs,
            args.w,
            tcc=args.tcc,
            mm=args.mm,
            filter=args.filter,
            threads=args.t,
            memory=args.m,
            overwrite=args.overwrite,
            loom=args.loom,
            h5ad=args.h5ad,
            cellranger=args.cellranger,
            report=args.report,
            inspect=not args.no_inspect,
            nucleus=args.workflow == 'nucleus',
            temp_dir=temp_dir
        )
    else:
        if args.workflow == 'kite:10xFB' and args.x.upper() != '10XV3':
            parser.error(
                '`kite:10xFB` workflow is only supported with technology `10XV3`'
            )

        # Smart-seq
        if args.x.upper() == 'SMARTSEQ':
            if args.dry_run:
                parser.error(f'Technology `{args.x}` does not support dry run.')

            if args.workflow != 'standard':
                parser.error(
                    f'TECHNOLOGY `{args.x}` only supports `standard` workflow.'
                )

            # Check for ignored arguments. (i.e. arguments either not supported or
            # not yet implemented)
            ignored = ['w', 'tcc', 'mm', 'filter', 'cellranger', 'report']
            for arg in ignored:
                if getattr(args, arg):
                    logger.warning(
                        f'Argument `{arg}` is not supported for technology `{args.x}`. This argument will be ignored.'
                    )

            # Check if batch TSV was provided.
            batch_path = None
            if len(args.fastqs) == 1:
                try:
                    with open_as_text(args.fastqs[0], 'r') as f:
                        if not f.readline().startswith('@'):
                            batch_path = args.fastqs[0]
                except Exception:
                    pass

            cells = {}
            if batch_path:
                with open(batch_path, 'r') as f:
                    for line in f:
                        if line.isspace() or line.startswith('#'):
                            continue
                        cell_id, fastq_1, fastq_2 = line.strip().split('\t')
                        if cell_id in cells:
                            parser.error(
                                f'Found duplicate cell ID {cell_id} in {batch_path}.'
                            )
                        cells[cell_id] = (fastq_1, fastq_2)
            else:
                # Allow glob notation for fastqs.
                fastqs = []
                for expr in args.fastqs:
                    fastqs.extend(glob.glob(expr))
                if len(fastqs) % 2:
                    logger.warning(f'{len(fastqs)} FASTQs found. ')
                fastqs = sorted(set(fastqs))
                for cell_id, i in enumerate(range(0, len(fastqs), 2)):
                    fastq_1, fastq_2 = fastqs[i], (
                        fastqs[i + 1] if i + 1 < len(fastqs) else ''
                    )
                    cells[cell_id] = (fastq_1, fastq_2)
            logger.info('Found the following FASTQs:')
            fastq_pairs = []
            cell_ids = []
            for cell_id, (fastq_1, fastq_2) in cells.items():
                logger.info(f'        {cell_id}    {fastq_1}  {fastq_2}')
                cell_ids.append(cell_id)
                fastq_pairs.append((fastq_1, fastq_2))
                if not fastq_1 or not fastq_2:
                    parser.error(
                        f'Single-end reads are currently not supported with technology `{args.x}`.'
                    )
            count_smartseq(
                args.i,
                args.g,
                args.x,
                args.o,
                fastq_pairs,
                cell_ids=cell_ids,
                threads=args.t,
                memory=args.m,
                overwrite=args.overwrite,
                loom=args.loom,
                h5ad=args.h5ad,
                temp_dir=temp_dir
            )
        else:
            count(
                args.i,
                args.g,
                args.x,
                args.o,
                args.fastqs,
                args.w,
                tcc=args.tcc,
                mm=args.mm,
                filter=args.filter,
                kite='kite' in args.workflow,
                FB='10xFB' in args.workflow,
                threads=args.t,
                memory=args.m,
                overwrite=args.overwrite,
                loom=args.loom,
                h5ad=args.h5ad,
                cellranger=args.cellranger,
                report=args.report,
                inspect=not args.no_inspect,
                temp_dir=temp_dir
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
    workflow = 'standard'
    for i, arg in enumerate(sys.argv):
        if arg.startswith('--workflow'):
            if '=' in arg:
                workflow = arg[arg.index('=') + 1:].strip('\'\"')
            else:
                workflow = sys.argv[i + 1]
            break

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
        help=(
            'Path to the kallisto index to be constructed. '
            'If `-n` is also specified, this is the prefix for the n indices to construct.'
        ),
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
        help=(
            '[Optional with -d] Path to the cDNA FASTA (lamanno, nucleus) '
            'or mismatch FASTA (kite) to be generated '
        ),
        type=str,
        required='-d' not in sys.argv
    )
    required_lamanno = parser_ref.add_argument_group(
        'required arguments for `lamanno` and `nucleus` workflows'
    )
    required_lamanno.add_argument(
        '-f2',
        metavar='FASTA',
        help='Path to the intron FASTA to be generated',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
    )
    required_lamanno.add_argument(
        '-c1',
        metavar='T2C',
        help='Path to generate cDNA transcripts-to-capture',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
    )
    required_lamanno.add_argument(
        '-c2',
        metavar='T2C',
        help='Path to generate intron transcripts-to-capture',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
    )

    parser_ref.add_argument(
        '-n',
        metavar='N',
        help=(
            'Number of files to split the index into. If this option is specified, '
            'the FASTA that is normally used to create an index is split into '
            '`N` approximately-equal parts. Each of these FASTAs are indexed separately. '
            'When using this option with `--workflow lamanno`, the intron FASTA '
            'is split into N-1 approximately-equal parts and indexed, while the '
            'cDNA FASTA is not split and indexed.'
        ),
        type=int,
        default=1,
        required=False
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
        '-k',
        metavar='K',
        help=(
            'Use this option to override the k-mer length of the index. '
            'Usually, the k-mer length automatically calculated by `kb` provides '
            'the best results.'
        ),
        type=int,
        default=None,
        required=False
    )
    parser_ref.add_argument(
        '--workflow',
        help=(
            'Type of workflow to prepare files for. '
            'Use `lamanno` for RNA velocity based on La Manno et al. 2018 logic. '
            'Use `nucleus` for RNA velocity on single-nucleus RNA-seq reads. '
            'Use `kite` for feature barcoding. (default: standard)'
        ),
        type=str,
        default='standard',
        choices=['standard', 'lamanno', 'nucleus', 'kite']
    )
    parser_ref.add_argument(
        '--overwrite',
        help='Overwrite existing kallisto index',
        action='store_true'
    )
    parser_ref.add_argument(
        'fasta',
        help='Genomic FASTA file(s), comma-delimited',
        type=str,
        nargs=None if '-d' not in sys.argv and workflow != 'kite' else '?'
    )
    parser_ref.add_argument(
        'gtf',
        help='Reference GTF file(s), comma-delimited',
        type=str,
        nargs=None if '-d' not in sys.argv and workflow != 'kite' else '?'
    )
    parser_ref.add_argument(
        'feature',
        help=(
            '[`kite` workflow only] Path to TSV containing barcodes and feature names.'
        ),
        type=str,
        nargs=None if '-d' not in sys.argv and workflow == 'kite' else '?'
    )

    # Hidden options.
    parser_ref.add_argument(
        '--no-mismatches', help=argparse.SUPPRESS, action='store_true'
    )
    parser_ref.add_argument('--flank', help=argparse.SUPPRESS, type=int)

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
    workflow = 'standard'
    for i, arg in enumerate(sys.argv):
        if arg.startswith('--workflow'):
            if '=' in arg:
                workflow = arg[arg.index('=') + 1:].strip('\'\"')
            else:
                workflow = sys.argv[i + 1]
            break

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
        help='Path to kallisto index/indices, comma-delimited',
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
        required=True,
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
            'a pre-packaged whitelist is used. Otherwise, or if \'None\', is '
            'provided, the bustools whitelist command is used. '
            '(`kb --list` to view whitelists)'
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
    parser_count.add_argument(
        '--workflow',
        help=(
            'Type of workflow. '
            'Use `lamanno` for RNA velocity based on La Manno et al. 2018 logic. '
            'Use `nucleus` for RNA velocity on single-nucleus RNA-seq reads. '
            'Use `kite` for feature barcoding. '
            'Use `kite:10xFB` for 10x Genomics Feature Barcoding technology. '
            '(default: standard)'
        ),
        type=str,
        default='standard',
        choices=['standard', 'lamanno', 'nucleus', 'kite', 'kite:10xFB']
    )

    count_group = parser_count.add_mutually_exclusive_group()
    count_group.add_argument(
        '--mm',
        help='Include reads that pseudoalign to multiple genes.',
        action='store_true'
    )
    count_group.add_argument(
        '--tcc',
        help='Generate a TCC matrix instead of a gene count matrix.',
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
    required_lamanno = parser_count.add_argument_group(
        'required arguments for `lamanno` and `nucleus` workflows'
    )
    required_lamanno.add_argument(
        '-c1',
        metavar='T2C',
        help='Path to cDNA transcripts-to-capture',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
        or any(arg in sys.argv for arg in {'--lamanno', '--nucleus'})
    )
    required_lamanno.add_argument(
        '-c2',
        metavar='T2C',
        help='Path to intron transcripts-to-captured',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
        or any(arg in sys.argv for arg in {'--lamanno', '--nucleus'})
    )
    parser_count.add_argument(
        '--overwrite',
        help='Overwrite existing output.bus file',
        action='store_true'
    )
    parser_count.add_argument('--dry-run', help='Dry run', action='store_true')

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
    parser_count.add_argument(
        '--cellranger',
        help='Convert count matrices to cellranger-compatible format',
        action='store_true'
    )
    report_group = parser_count.add_mutually_exclusive_group()
    report_group.add_argument(
        '--report',
        help=(
            'Generate a HTML report containing run statistics and basic plots. '
            'Using this option may cause kb to use more memory than specified '
            'with the `-m` option. It may also cause it to crash due to memory.'
        ),
        action='store_true'
    )
    report_group.add_argument(
        '--no-inspect', help=argparse.SUPPRESS, action='store_true'
    )
    parser_count.add_argument(
        '--no-validate', help=argparse.SUPPRESS, action='store_true'
    )
    parser_count.add_argument(
        'fastqs',
        help=(
            'FASTQ files. For technology `SMARTSEQ`, all input FASTQs are '
            'alphabetically sorted by path and paired in order, and cell IDs '
            'are assigned as incrementing integers starting from zero. A single '
            'batch TSV with cell ID, read 1, and read 2 as columns can be '
            'provided to override this behavior.'
        ),
        nargs='+'
    )
    return parser_count


@logger.namespaced('main')
def main():
    """Command-line entrypoint.
    """
    # Get prepackaged kallisto and bustools paths.
    kallisto_path = get_kallisto_binary_path()
    bustools_path = get_bustools_binary_path()

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
        '--tmp',
        metavar='TMP',
        help='Override default temporary directory',
        type=str
    )
    parent.add_argument(
        '--keep-tmp',
        help='Do not delete the tmp directory',
        action='store_true'
    )
    parent.add_argument(
        '--verbose', help='Print debugging information', action='store_true'
    )
    parent.add_argument(
        '--kallisto',
        help=f'Path to kallisto binary to use (default: {kallisto_path})',
        type=str,
        default=kallisto_path
    )
    parent.add_argument(
        '--bustools',
        help=f'Path to bustools binary to use (default: {bustools_path})',
        type=str,
        default=bustools_path
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

    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    # Validation
    if 'no_validate' in args and args.no_validate:
        logger.warning((
            'File validation is turned off. '
            'This may lead to corrupt/empty output files.'
        ))
        no_validate()

    if 'dry_run' in args:
        # Dry run can not be specified with matrix conversion.
        if args.dry_run and (args.loom or args.h5ad):
            raise parser.error(
                '--dry-run can not be used with --loom or --h5ad'
            )

        if args.dry_run:
            logging.disable(level=logging.CRITICAL)
            set_dry()

    logger.debug('Printing verbose output')

    # Set binary paths
    set_kallisto_binary_path(args.kallisto)
    set_bustools_binary_path(args.bustools)

    logger.debug(f'kallisto binary located at {get_kallisto_binary_path()}')
    logger.debug(f'bustools binary located at {get_bustools_binary_path()}')

    temp_dir = args.tmp or (
        os.path.join(args.o, TEMP_DIR) if 'o' in args else TEMP_DIR
    )
    # Check if temp_dir exists and exit if it does.
    # This is so that kb doesn't accidently use an existing directory and
    # delete it afterwards.
    if os.path.exists(temp_dir):
        parser.error(
            f'Temporary directory `{temp_dir}` exists! Is another instance running? '
            'Either remove the existing directory or use `--tmp` to specify a '
            'different temporary directory.'
        )

    logger.debug(f'Creating `{temp_dir}` directory')
    make_directory(temp_dir)
    try:
        logger.debug(args)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            COMMAND_TO_FUNCTION[args.command](parser, args, temp_dir=temp_dir)
    except Exception:
        if is_dry():
            raise
        logger.exception('An exception occurred')
    finally:
        # Always clean temp dir
        if not args.keep_tmp:
            logger.debug(f'Removing `{temp_dir}` directory')
            remove_directory(temp_dir)
