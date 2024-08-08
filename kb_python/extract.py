import os
import shutil
from typing import Dict, List, Optional, Union

import pandas as pd

from typing_extensions import Literal

from .utils import make_directory, get_bustools_binary_path, run_executable

from .count import kallisto_bus, bustools_capture, bustools_sort, stream_fastqs

from .logging import logger


def bustools_extract(
    sorted_bus_path: str,
    out_path: str,
    fastqs: Union[str, List[str]],
    num_fastqs: int,
) -> Dict[str, str]:
    """
    Extract reads from a BUS file using bustools.
    Args:
        sorted_bus_path: path to a BUS file sorted by flag using bustools sort --flag
        out_path: Output directory for FASTQ files
        fastqs: FASTQ file(s) from which to extract reads
        num_fastqs: Number of FASTQ file(s) per run

    Returns:
        Dictionary containing path to generated BUS file
    """
    logger.info(f"Extracting BUS file {sorted_bus_path} to {out_path}")
    command = [get_bustools_binary_path(), 'extract']

    if not isinstance(fastqs, list):
        fastqs = [fastqs]

    command += ['-o', out_path]
    command += ['-f', ",".join(fastqs)]
    command += ['-N', num_fastqs]
    command += [sorted_bus_path]
    run_executable(command)
    return {"bus": out_path}


@logger.namespaced('extract')
def extract(
    fastq: str,
    index_path: str,
    targets: list[str],
    out_dir: str,
    target_type: Literal['gene', 'transcript'],
    extract_all: False,
    t2g_path: Optional[str] = None,
    temp_dir: str = 'tmp',
    threads: int = 8,
    aa: bool = False,
    strand: Optional[Literal['unstranded', 'forward', 'reverse']] = None,
    numreads: Optional[int] = None,
):
    """
    Extracts sequencing reads that were pseudo-aligned to an index for specific genes/transcripts.
    Note: Multimapped reads will also be extracted.
    
    fastq: Single fastq file containing sequencing reads
    index_path: Path to kallisto index
    targets: Gene or transcript names for which to extract the raw reads that align to the index
    out_dir: Path to output directory
    target_type: 'gene' (default) or 'transcript' -> Defines whether targets are gene or transcript names
    extract_all: Extracts reads for all genes, defaults to `False`. Set targets = None when using extract_all.
    t2g_path: Path to transcript-to-gene mapping file (required when target_type = gene or extract_all = True)
    temp_dir: Path to temporary directory, defaults to `tmp`
    threads: Number of threads to use, defaults to `8`
    aa: Align to index generated from a FASTA-file containing amino acid sequences, defaults to `False`
    strand: Strandedness, defaults to `None`
    numreads: Maximum number of reads to process from supplied input

    Returns:
    Raw reads that were pseudo-aligned to the index by kallisto for each specified gene/transcript.
    """
    if target_type not in ["gene", "transcript"]:
        raise ValueError(
            f"target_type must be 'gene' or 'transcript', not {target_type}"
        )

    if (target_type == "gene" or extract_all) and (t2g_path is None):
        raise ValueError(
            "t2g_path must be provided if target_type is 'gene' or extract_all is True"
        )

    make_directory(out_dir)

    fastq = stream_fastqs([fastq], temp_dir=temp_dir)

    logger.info("Performing alignment using kallisto...")

    kallisto_bus(
        fastqs=fastq,
        index_path=index_path,
        technology="bulk",
        out_dir=temp_dir,
        threads=threads,
        n=True,
        paired=False,
        aa=aa,
        strand=strand,
        numreads=numreads
    )

    logger.info("Alignment complete. Beginning extraction of reads using bustools...")
    if target_type == "gene" or extract_all:
        # Read t2g to find all transcripts associated with a gene/mutant ID
        with open(t2g_path, "r") as t2g_file:
            lines = t2g_file.readlines()
        t2g_df = pd.DataFrame()
        t2g_df["transcript"] = [line.split("\t")[0] for line in lines]
        t2g_df["gene_id"] = [
            line.split("\t")[1].replace("\n", "") for line in lines
        ]

        # If extract_all is True, extract pseudo-aligned reads for all genes
        if extract_all:
            targets = t2g_df["gene_id"].values

        g2ts = {
            gid: t2g_df[t2g_df["gene_id"] == gid]["transcript"].values.tolist()
            for gid in targets
        }

    ecmap = os.path.join(temp_dir, "matrix.ec")
    txnames = os.path.join(temp_dir, "transcripts.txt")
    bus_in = os.path.join(temp_dir, "output.bus")

    for gid in targets:
        if target_type == "gene" or extract_all:
            transcripts = g2ts[gid]
        else:
            # if target_type==transcript, each transcript will be extracted individually
            transcripts = [gid]

        # Create temp txt file with transcript IDs to extract
        transcript_names_file = os.path.join(
            temp_dir, "pull_out_reads_transcript_ids_temp.txt"
        )
        with open(transcript_names_file, "w") as f:
            f.write("\n".join(transcripts))

        if target_type == "gene":
            logger.info(
                f"Extracting reads for following transcripts for gene ID {gid}: "
                + ", ".join(transcripts)
            )
        else:
            logger.info(
                f"Extracting reads for the following transcript: {gid}"
            )

        bus_out = os.path.join(temp_dir, f"output_extracted_{gid}.bus")
        bus_out_sorted = os.path.join(
            temp_dir, f"output_extracted_{gid}_sorted.bus"
        )

        try:
            # Capture records for this transcript ID
            bustools_capture(
                bus_path=bus_in,
                capture_path=transcript_names_file,
                ecmap_path=ecmap,
                txnames_path=txnames,
                capture_type="transcripts",
                out_path=bus_out,
                complement=False
            )
    
            # Extract records for this transcript ID from fastq
            # bustools_sort(bus_path=bus_out, flags=True, out_path=bus_out_sorted)

            # Omit sorting because 'bustools sort' has trouble when a flag column is present
            shutil.copyfile(bus_out, bus_out_sorted)
    
            extract_out_folder = os.path.join(out_dir, gid)
            bustools_extract(
                sorted_bus_path=bus_out_sorted,
                out_path=extract_out_folder,
                fastqs=fastq,
                num_fastqs=1
            )

        except Exception as e:
            logger.error(f"Extraction of reads unsuccessful for {gid} due to the following error:\n{e}")
