from urllib.parse import urlparse


def stream_batch(batch_path, temp_dir='tmp'):
    """Dry version of `count.stream_batch`.
    """
    with open(batch_path, 'r') as f:
        for line in f:
            if line.isspace() or line.startswith('#'):
                continue

            split = line.split('\t')
            if any(urlparse(fastq).scheme in ('http', 'https', 'ftp', 'ftps')
                   for fastq in split[1:]):
                raise Exception(
                    'Streaming remote FASTQs from a batch file is not dryable.'
                )
