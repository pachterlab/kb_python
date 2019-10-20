import logging
import os
import shutil

from .constants import COMBINED_FILENAME
from .utils import (
    concatenate_files,
    create_transcript_list,
    get_transcripts_from_fasta,
    run_executable,
)

logger = logging.getLogger(__name__)


def kallisto_index(fasta_path, index_path, k=31):
    command = ['kallisto', 'index', '-i', index_path, '-k', k, fasta_path]
    run_executable(command)
    return {'index': index_path}


def create_t2g(
        gtf_path, t2g_path, use_name=True, use_version=True, transcripts=None
):
    r = create_transcript_list(
        gtf_path,
        use_name=use_name,
        use_version=use_version,
        transcripts=None,
    )
    with open(t2g_path, 'w') as f:
        for tid in r:
            if use_name:
                f.write('{}\t{}\t{}\n'.format(tid, r[tid][0], r[tid][1]))
            else:
                f.write('{}\t{}\n'.format(tid, r[tid][0]))
    return {'t2g': t2g_path}


def create_t2c(fasta_path, t2c_path):
    transcripts = get_transcripts_from_fasta(fasta_path)
    with open(t2c_path, 'w') as f:
        f.write('\n'.join(transcripts))
    return {'t2c': t2c_path, 'transcripts': transcripts}


def ref(fasta_path, gtf_path, index_path, t2g_path, overwrite=False):
    results = {}
    t2g_result = create_t2g(gtf_path, t2g_path)
    results.update(t2g_result)
    if not os.path.exists(index_path) or overwrite:
        index_result = kallisto_index(fasta_path, index_path)
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )
    return results


def ref_velocity(
        cdna_path,
        intron_path,
        gtf_path,
        index_path,
        t2g_path,
        cdna_t2c_path,
        intron_t2c_path,
        temp_dir='tmp',
        keep_temp=False,
        overwrite=False
):
    # Make temporary directory.
    os.makedirs(temp_dir, exist_ok=True)
    results = {}

    cdna_t2c_result = create_t2c(cdna_path, cdna_t2c_path)
    results.update({'cdna_t2c': cdna_t2c_result['t2c']})
    intron_t2c_result = create_t2c(intron_path, intron_t2c_path)
    intron_without_prefix = [
        t for t in intron_t2c_result['transcripts'] if not t.endswith('-I')
    ]
    if intron_without_prefix:
        logger.warning((
            'Found {} intron transcript IDs that do not have the "-I" prefix. '
            'Downstream analyses may be compromised.'
        ).format(len(intron_without_prefix)))
    results.update({'intron_t2c': intron_t2c_result['t2c']})

    transcripts = set(
        cdna_t2c_result['transcripts'] + intron_t2c_result['transcripts']
    )
    t2g_result = create_t2g(gtf_path, t2g_path, transcripts=transcripts)
    results.update(t2g_result)

    if not os.path.exists(index_path) or overwrite:
        os.makedirs(temp_dir, exist_ok=True)
        combined_path = concatenate_files(
            cdna_path,
            intron_path,
            out_path=os.path.join(temp_dir, COMBINED_FILENAME),
            temp_dir=temp_dir
        )
        index_result = kallisto_index(combined_path, index_path)
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )

    # Remove temporary directory.
    if not keep_temp:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return results
