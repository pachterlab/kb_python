from collections import namedtuple

# Technology to file position mapping
Technology = namedtuple('Technology', ['name', 'nfiles', 'seq_pos', 'umi_pos'])
TECHNOLOGIES = [
    Technology('10XV1', 3, 0, 1),
    Technology('10XV2', 2, 1, 0),
    Technology('10XV3', 2, 1, 0),
    Technology('CELSEQ', 2, 1, 0),
    Technology('CELSEQ2', 2, 1, 0),
    Technology('DROPSEQ', 2, 1, 0),
    Technology('INDROPS', 2, 1, 0),
    Technology('SCRUBSEQ', 2, 1, 0),
    Technology('SURECELL', 2, 1, 0),
]
