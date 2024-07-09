import os
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
    index_path: str,
    technology: str,
    out_dir: str,
    fastq: str,
    targets: list[str],
    target_type: Literal['gene', 'transcript'],
    t2g_path: Optional[str] = None,
    temp_dir: str = 'tmp',
    threads: int = 8,
    aa: bool = False,
    strand: Optional[Literal['unstranded', 'forward', 'reverse']] = None,
    numreads: Optional[int] = None,
):
    if target_type not in ["gene", "transcript"]:
        raise ValueError(
            f"target_type must be 'gene' or 'transcript', not {target_type}"
        )

    if (target_type == "gene") == (t2g_path is None):
        raise ValueError(
            "t2g_path must be provided if and only if target_type is 'gene'"
        )

    make_directory(out_dir)

    fastq = stream_fastqs([fastq], temp_dir=temp_dir)[0]

    kallisto_bus(
        fastqs=fastq,
        index_path=index_path,
        technology=technology,
        out_dir=out_dir,
        threads=threads,
        n=True,
        paired=False,
        aa=aa,
        strand=strand,
        numreads=numreads
    )

    if target_type == "gene":
        # Extract reads using bustools
        # Read t2g to find all transcripts associated with a gene/mutant ID
        with open(t2g_path, "r") as t2g_file:
            lines = t2g_file.readlines()
        t2g_df = pd.DataFrame()
        t2g_df["transcript"] = [line.split("\t")[0] for line in lines]
        t2g_df["gene_id"] = [
            line.split("\t")[1].replace("\n", "") for line in lines
        ]
        gene_ids = targets
        g2ts = {
            gid: t2g_df[t2g_df["gene_id"] == gid]["transcript"].values.tolist()
            for gid in gene_ids
        }
    else:
        gene_ids = ["transcript"]
        g2ts = {gene_ids[0]: targets}

    ecmap = os.path.join(out_dir, "matrix.ec")
    txnames = os.path.join(out_dir, "transcripts.txt")
    bus_in = os.path.join(out_dir, "output.bus")

    for gid in gene_ids:
        transcripts = g2ts[gid]

        # Create temp txt file with transcript IDs to extract
        transcript_names_file = os.path.join(
            temp_dir, "pull_out_reads_transcript_ids_temp.txt"
        )
        with open(transcript_names_file, "w") as f:
            f.write("\n".join(transcripts))

        if target_type == "gene":
            logger.info(
                f"Extracting reads from the following transcripts for gene ID {gid}: "
                + ", ".join(transcripts)
            )
        else:
            logger.info(
                "Extracting reads from the following transcripts: " +
                ", ".join(transcripts)
            )

        bus_out = os.path.join(out_dir, f"output_extracted_{gid}.bus")
        bus_out_sorted = os.path.join(
            out_dir, f"output_extracted_{gid}_sorted.bus"
        )

        # Capture records for this transcript ID
        bustools_capture(
            bus_path=bus_in,
            capture_path=transcript_names_file,
            ecmap_path=ecmap,
            txnames_path=txnames,
            capture_type="transcripts",
            out_path=bus_out
        )

        # Extract records for this transcript ID from fastq
        bustools_sort(bus_path=bus_out, flags=True, out_path=bus_out_sorted)

        extract_out_folder = os.path.join(out_dir, gid)
        bustools_extract(
            sorted_bus_path=bus_out_sorted,
            out_path=extract_out_folder,
            fastqs=fastq,
            num_fastqs=1
        )
