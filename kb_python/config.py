import os
from collections import namedtuple

# Technology to file position mapping
Technology = namedtuple('Technology', ['name', 'nfiles', 'seq_pos', 'umi_pos', 'whitelist_archive', 'whitelist_filename'])
WHITELIST_DIR = 'whitelists'
TECHNOLOGIES = [
    Technology('10XV1', 3, 0, 1, '10xv1_whitelist.tar.gz', '10xv1_whitelist.txt'),
    Technology('10XV2', 2, 1, 0, '10xv2_whitelist.tar.gz', '10xv2_whitelist.txt'),
    Technology('10XV3', 2, 1, 0, '10xv3_whitelist.tar.gz', '10xv3_whitelist.txt'),
    Technology('CELSEQ', 2, 1, 0, None, None),
    Technology('CELSEQ2', 2, 1, 0, None, None),
    Technology('DROPSEQ', 2, 1, 0, None, None),
    Technology('INDROPS', 2, 1, 0, None, None),
    Technology('SCRUBSEQ', 2, 1, 0, None, None),
    Technology('SURECELL', 2, 1, 0, None, None),
]
TECHNOLOGIES_MAPPING = {t.name: t for t in TECHNOLOGIES}
