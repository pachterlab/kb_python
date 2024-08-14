import os
from typing import Dict, List, Optional, Union

from Bio import SeqIO
import pandas as pd
import numpy as np

from typing_extensions import Literal

from .utils import make_directory, get_bustools_binary_path, run_executable

from .count import kallisto_bus, bustools_capture, bustools_sort, stream_fastqs

from .logging import logger

def bustools_text(
    bus_path: str,
    out_path: str,
    flags: bool = False,
):
    """Runs `bustools text`.

    Args:
        bus_path: Path to BUS file to convert to text format
        out_dir: Path to output txt file path
        flags: Whether to include the flags columns
    """
    # logger.info('Converting BUS file {} to {}'.format(bus_path, out_path))
    command = [get_bustools_binary_path(), 'text']
    command += ['-o', out_path]
    if flags:
        command += ['--flags']
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}

def bustools_fromtext(
    txt_path: str,
    out_path: str
):
    """Runs `bustools fromtext`.

    Args:
        bus_path: Path to text file to convert to BUS format
        out_dir: Path to output BUS file
    """
    # logger.info('Creating BUS file {} from {}'.format(out_path, bus_path))
    command = [get_bustools_binary_path(), 'fromtext']
    command += ['-o', out_path]
    command += [txt_path]
    run_executable(command)
    return {'bus': out_path}

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

def read_headers_from_fastq(fastq_file):
    """
    Reads headers from a FASTQ file and returns a set of headers.
    """
    headers = set()
    with open(fastq_file, 'r') as file:
        for record in SeqIO.parse(file, 'fastq'):
            headers.add(record.id)
    return headers

def extract_matching_reads_by_header(input_fastq, reference_fastq, output_fastq):
    """
    Extracts reads from the input FASTQ file that are present in the reference FASTQ file
    based on headers and writes them to the output FASTQ file.
    """
    # Read headers from the reference FASTQ file
    reference_headers = read_headers_from_fastq(reference_fastq)
    
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        # Create a SeqIO writer for the output FASTQ file
        writer = SeqIO.write(
            (record for record in SeqIO.parse(infile, 'fastq') if record.id in reference_headers),
            outfile,
            'fastq'
        )
        logger.info(f'Number of unmapped reads written to {output_fastq}: {writer}')

@logger.namespaced('extract')
def extract(
    fastq: str,
    index_path: str,
    targets: list[str],
    out_dir: str,
    target_type: Literal['gene', 'transcript'],
    extract_all: bool = False,
    extract_all_fast: bool = False,
    extract_all_unmapped: bool = False,
    mm: bool = False,
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
    extract_all: Extracts reads for all genes or transcripts (as defined in target_type), defaults to `False`. Might take a long time to run when the reference index contains a large number of genes. Set targets = None when using extract_all.
    extract_all_fast: Extracts all pseudo-aligned reads, defaults to `False`. Does not break down output by gene/transcript. Set targets = None when using extract_all_fast.
    extract_all_unmapped: Extracts all unmapped reads, defaults to `False`. Set targets = None when using extract_all_unmapped.
    mm: Also extract reads that multi-mapped to several genes, defaults to `False`
    t2g_path: Path to transcript-to-gene mapping file (required when target_type = gene or extract_all = True)
    temp_dir: Path to temporary directory, defaults to `tmp`
    threads: Number of threads to use, defaults to `8`
    aa: Align to index generated from a FASTA-file containing amino acid sequences, defaults to `False`
    strand: Strandedness, defaults to `None`
    numreads: Maximum number of reads to process from supplied input

    Returns:
    Raw reads that were pseudo-aligned to the index by kallisto for each specified gene/transcript.
    """
    if sum([extract_all, extract_all_fast, extract_all_unmapped]) > 1:
        raise ValueError(
            f"extract_all, extract_all_fast, and/or extract_all_unmapped cannot be used simultaneously"
        )

    if targets is None and not (extract_all or extract_all_fast or extract_all_unmapped):
        raise ValueError(
            f"targets must be provided (unless extract_all, extract_all_fast, or extract_all_unmapped are used to extract all reads)"
        )

    if targets and (extract_all or extract_all_fast or extract_all_unmapped):
        logger.warning(
            f"targets will be ignored since extract_all, extract_all_fast, or extract_all_unmapped is activated which will extract all reads"
        )

    if target_type not in ["gene", "transcript"]:
        raise ValueError(
            f"target_type must be 'gene' or 'transcript', not {target_type}"
        )

    if (not mm or (target_type == "gene" and not (extract_all_fast or extract_all_unmapped)) or extract_all) and (t2g_path is None):
        raise ValueError(
            "t2g_path must be provided if mm flag is not provided, target_type is 'gene' (and extract_all_fast and extract_all_unmapped are False), OR extract_all is True"
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

    txnames = os.path.join(temp_dir, "transcripts.txt")
    bus_in = os.path.join(temp_dir, "output.bus")

    if mm:
        ecmap = os.path.join(temp_dir, "matrix.ec")

    # Remove multimapped reads from matrix and bus files
    else:
        logger.info("Removing equivalence classes with multi-mapped reads from the matrix.ec and BUS files")

        ecmap_mm = os.path.join(temp_dir, "matrix.ec")
        ecmap = os.path.join(temp_dir, "matrix_no_mm.ec")

        # Read t2g to find all transcripts associated with a gene/mutant ID
        with open(t2g_path, "r") as t2g_file:
            lines = t2g_file.readlines()
        t2g_df = pd.DataFrame()
        t2g_df["transcript"] = [line.split("\t")[0] for line in lines]
        t2g_df["gene_id"] = [
            line.split("\t")[1].replace("\n", "") for line in lines
        ]

        with open(txnames) as f:
            txs = f.read().splitlines()

        ec_df = pd.read_csv(ecmap_mm, sep="\t", header=None)
        # List to save multimapped ecs
        ecs_mm = []
        for index, row in ec_df.iterrows():
            # Get transcript IDs that mapped to this ec
            mapped_txs = np.array(txs)[np.array(row[1].split(",")).astype(int)]

            # Check if transcript IDs belong to one or more genes
            if len(set(t2g_df[t2g_df["transcript"].isin(mapped_txs)]["gene_id"].values)) > 1:
                ecs_mm.append(row[0])

        # Write new matrix.ec file excluding mm ecs
        new_ec = ec_df[~ec_df[0].isin(ecs_mm)]
        new_ec.to_csv(ecmap, sep='\t', index=False, header=None)

        logger.debug("Finished removing equivalence classes with multimapped reads from matrix.ec")

        ## Remove mm ecs from bus file
        bus_txt = os.path.join(temp_dir, "output.bus.txt")
        bus_txt_no_mm = os.path.join(temp_dir, "output_no_mm.bus.txt")
        bus_no_mm = os.path.join(temp_dir, "output_no_mm.bus")

        # Convert bus to txt file
        bustools_text(
            bus_path = bus_in,
            out_path = bus_txt,
            flags = True
        )

        # Remove mm ecs
        bus_df = pd.read_csv(bus_txt, sep="\t", header=None)
        new_bus_df = bus_df[~bus_df[2].isin(ecs_mm)]
        new_bus_df.to_csv(bus_txt_no_mm, sep='\t', index=False, header=None)

        # Convert back to bus format
        bustools_fromtext(
            txt_path = bus_txt_no_mm,
            out_path = bus_no_mm            
        )

        bus_in = bus_no_mm

        logger.debug("Finished removing equivalence classes with multimapped reads from BUS file")

    if extract_all_fast or extract_all_unmapped:
        bus_out_sorted = os.path.join(
            temp_dir, "output_extracted_sorted.bus"
        )

        try:
            # Extract records for this transcript ID from fastq
            bustools_sort(bus_path=bus_in, flags=True, out_path=bus_out_sorted)
    
            extract_out_folder = os.path.join(out_dir, "all")
            bustools_extract(
                sorted_bus_path=bus_out_sorted,
                out_path=extract_out_folder,
                fastqs=fastq,
                num_fastqs=1
            )

        except Exception as e:
            logger.error(f"Extraction of reads unsuccessful due to the following error:\n{e}")

    if extract_all_unmapped:
        # Save unmapped reads in a separate fastq file
        unmapped_fastq = os.path.join(out_dir, "all_unmapped")
        extract_matching_reads_by_header(extract_out_folder, fastq, unmapped_fastq)

    else:
        if target_type == "gene" or extract_all:
            # Read t2g to find all transcripts associated with a gene/mutant ID
            with open(t2g_path, "r") as t2g_file:
                lines = t2g_file.readlines()
            t2g_df = pd.DataFrame()
            t2g_df["transcript"] = [line.split("\t")[0] for line in lines]
            t2g_df["gene_id"] = [
                line.split("\t")[1].replace("\n", "") for line in lines
            ]
    
            if extract_all:
                if target_type == "gene":
                    # Set targets to all genes
                    targets = list(set(t2g_df["gene_id"].values))
                    g2ts = {
                            gid: t2g_df[t2g_df["gene_id"] == gid]["transcript"].values.tolist()
                            for gid in targets
                        }
    
                else:
                    # Set targets to all transcripts
                    targets = list(set(t2g_df["transcript"].values))
    
            else:
                g2ts = {
                        gid: t2g_df[t2g_df["gene_id"] == gid]["transcript"].values.tolist()
                        for gid in targets
                    }
    
        for gid in targets:
            if target_type == "gene":
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
                # Capture records for this gene
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
                bustools_sort(bus_path=bus_out, flags=True, out_path=bus_out_sorted)
        
                extract_out_folder = os.path.join(out_dir, gid)
                bustools_extract(
                    sorted_bus_path=bus_out_sorted,
                    out_path=extract_out_folder,
                    fastqs=fastq,
                    num_fastqs=1
                )
    
            except Exception as e:
                logger.error(f"Extraction of reads unsuccessful for {gid} due to the following error:\n{e}")
