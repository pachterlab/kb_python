# Implementation Summary: GitHub Issues #295 and #296

## Overview
This implementation addresses two feature requests for kb_python:
- **Issue #295**: Add output mode compatible with Read10X() in Seurat (spliced/unspliced separated as in CellRanger format)
- **Issue #296**: Enable gzip compression of output files and optional cleanup of BUS files to reduce disk usage

## Changes Made

### 1. New Command-Line Flags (main.py)

Added three new flags to the `kb count` command:

#### `--cellranger`
- **Purpose**: Convert count matrices to CellRanger-compatible format
- **Auto-enabled features**:
  - **Gzip compression**: Automatically compresses output matrices (can be disabled with `--no-gzip`)
  - **CellRanger-style directories**: For nac/lamanno workflows, automatically creates `spliced/` and `unspliced/` subdirectories
- **Files created**: `matrix.mtx.gz`, `barcodes.tsv.gz`, `genes.tsv.gz` (or `features.tsv.gz`)
- **Usage**: `kb count --cellranger ...` (gzip is automatic)
- **Disable gzip**: `kb count --cellranger --no-gzip ...`

#### `--gzip`
- **Purpose**: Manually enable gzip compression (not needed if using `--cellranger`)
- **Default**: False (but auto-enabled with `--cellranger`)
- **Usage**: `kb count --gzip ...` (without --cellranger)

#### `--no-gzip`
- **Purpose**: Disable automatic gzip compression when using `--cellranger`
- **Default**: False
- **Usage**: `kb count --cellranger --no-gzip ...`

#### `--delete-bus`
- **Purpose**: Delete intermediate BUS files after successful count to save disk space
- **Files deleted**: All `.bus` files generated during processing
- **Default**: False
- **Usage**: `kb count --delete-bus ...`
- **Safety**: Only deletes files after successful completion of count operation

### CellRanger-Style Output Structure

When using `--cellranger` with nac/lamanno workflows, the output is automatically organized as:
```
counts_unfiltered/
  spliced/
    matrix.mtx.gz
    barcodes.tsv.gz
    genes.tsv.gz
  unspliced/
    matrix.mtx.gz
    barcodes.tsv.gz
    genes.tsv.gz
  cellranger_ambiguous/
    matrix.mtx.gz
    barcodes.tsv.gz
    genes.tsv.gz
```

This structure is compatible with Seurat's `Read10X()` function and works with both filtered and unfiltered outputs.

### 2. Function Updates

#### `matrix_to_cellranger()` (count.py)
- **New parameter**: `gzip: bool = False`
- **Behavior**: 
  - When `gzip=True`, outputs `.gz` compressed versions of all files
  - Uses `compress_gzip()` from utils for matrix file compression
  - Uses `open_as_text()` for text files (barcodes, genes) to handle gzip automatically
- **Backward compatible**: Default behavior unchanged

#### `count()` function (count.py)
- **New parameters**: 
  - `gzip: bool = False`
  - `delete_bus: bool = False`
- **Changes**:
  - Passes `gzip` parameter to all `matrix_to_cellranger()` calls
  - Passes `gzip` parameter to `filter_with_bustools()`
  - Implements BUS file cleanup at end of function when `delete_bus=True`
  - Collects all BUS file paths from results and deletes them

#### `count_nac()` function (count.py)
- **New parameters**: 
  - `gzip: bool = False`
  - `cellranger_style: bool = False`
  - `delete_bus: bool = False`
- **Changes**:
  - Passes `gzip` parameter to all `matrix_to_cellranger()` calls
  - Implements `cellranger_style` logic:
    - For `processed` matrices (index 0): creates `spliced/` subdirectory
    - For `unprocessed` matrices (index 1): creates `unspliced/` subdirectory
    - For `ambiguous` matrices: uses default naming
  - Works for both filtered and unfiltered outputs
  - Implements BUS file cleanup at end of function when `delete_bus=True`

#### `filter_with_bustools()` function (count.py)
- **New parameter**: `gzip: bool = False`
- **Changes**: Passes `gzip` parameter to `matrix_to_cellranger()` call

### 3. Integration in main.py

#### `parse_count()` function
- Added logic to auto-enable gzip when `--cellranger` is used (unless `--no-gzip` is specified)
- Auto-enables cellranger-style directory structure for nac/lamanno workflows when `--cellranger` is used
- Passes appropriate parameters to count functions:
  - `gzip = (args.cellranger and not args.no_gzip) or args.gzip`
  - `cellranger_style = args.cellranger` (for count_nac)
  - `delete_bus = args.delete_bus`

## Usage Examples

### Example 1: Standard workflow with automatic gzip compression
```bash
kb count -i index.idx -g t2g.txt -x 10XV3 -o output/ \
  --cellranger \
  sample_R1.fastq.gz sample_R2.fastq.gz
```
Note: Gzip compression is now automatic! No need for `--gzip` flag.

### Example 2: NAC workflow with automatic CellRanger-style output
```bash
kb count -i index.idx -g t2g.txt -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  -x 10XV3 -o output/ --workflow=nac \
  --cellranger \
  sample_R1.fastq.gz sample_R2.fastq.gz
```
This automatically creates:
- `output/counts_unfiltered/spliced/` with spliced matrices (gzipped)
- `output/counts_unfiltered/unspliced/` with unspliced matrices (gzipped)

### Example 3: Disable gzip if needed
```bash
kb count -i index.idx -g t2g.txt -x 10XV3 -o output/ \
  --cellranger --no-gzip \
  sample_R1.fastq.gz sample_R2.fastq.gz
```

### Example 4: All features combined
```bash
kb count -i index.idx -g t2g.txt -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  -x 10XV3 -o output/ --workflow=nac \
  --cellranger --delete-bus \
  sample_R1.fastq.gz sample_R2.fastq.gz
```
This automatically enables gzip and creates CellRanger-style directories!

## Benefits

1. **Disk Space Savings**: 
   - Gzip compression reduces matrix file sizes by 70-90%
   - Automatically enabled with `--cellranger` flag
   - BUS file deletion frees up significant temporary storage
   - Particularly beneficial for large-scale processing pipelines

2. **Seurat Compatibility**: 
   - `--cellranger` automatically creates Seurat-compatible output
   - CellRanger-style directories created automatically for nac/lamanno workflows
   - Spliced/unspliced separation allows easy loading of separate assays
   - Compatible with existing R-based analysis workflows
   - No additional flags needed!

3. **Standard Compliance**: 
   - Gzipped outputs align with CellRanger format standards (enabled by default)
   - Scanpy and other tools natively support `.gz` files
   - No manual compression/decompression needed

4. **Simplified Usage**:
   - Just use `--cellranger` and everything is configured automatically
   - Gzip compression: automatic
   - CellRanger-style directories (for nac/lamanno): automatic
   - Can disable gzip with `--no-gzip` if needed

## Backward Compatibility

All changes are backward compatible:
- New flags are optional (default: False)
- Existing commands continue to work without modification
- Default behavior unchanged when new flags not specified
- Function signatures extended with optional parameters (defaults maintain old behavior)

## Testing Recommendations

1. Test gzip compression:
   - Verify `.gz` files are created
   - Verify compressed files can be read by Seurat/Scanpy
   - Compare file sizes before/after compression

2. Test CellRanger-style directories:
   - Verify `spliced/` and `unspliced/` subdirectories are created
   - Verify matrices contain correct data
   - Test with Seurat's `Read10X()` function

3. Test BUS file deletion:
   - Verify BUS files are deleted after successful run
   - Verify final matrices are still correct
   - Test that deletion doesn't occur if count fails

4. Test combinations:
   - All flags together
   - Each flag independently
   - With filtered and unfiltered outputs

## Files Modified

1. `/Users/benjamin/kb_python/kb_python/main.py`
   - Added `--cellranger`, `--gzip`, `--no-gzip`, and `--delete-bus` command-line arguments
   - Updated `parse_count()` to auto-enable gzip and cellranger-style when `--cellranger` is used
   - Logic: `use_gzip = (args.cellranger and not args.no_gzip) or args.gzip`
   - Logic: `use_cellranger_style = args.cellranger` (for nac/lamanno)

2. `/Users/benjamin/kb_python/kb_python/count.py`
   - Updated `matrix_to_cellranger()` with gzip support
   - Updated `count()` with gzip and delete_bus support
   - Updated `count_nac()` with gzip, cellranger_style, and delete_bus support
   - Updated `filter_with_bustools()` with gzip support
   - Added BUS file cleanup logic to count functions

## Notes

- **Simplified workflow**: `--cellranger` now automatically enables both gzip and cellranger-style directories
- The cellranger-style flag has been removed as it's redundant (automatic with `--cellranger`)
- Gzip compression uses existing `compress_gzip()` utility from ngs_tools
- BUS file deletion only occurs after successful completion of the count operation
- For smartseq3 technology, handles multiple BUS file suffixes correctly
- Use `--no-gzip` with `--cellranger` if you need uncompressed output
