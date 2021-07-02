# eCLIP-Peak

Pipeline for using IDR to identify a set of reproducible peaks given eClIP dataset with two or three replicates.

## Installation
- For Van Nostrand Lab

    The pipeline has already been installed. Activate its environment
    by issue the following command: 
    `source /storage/vannostrand/software/eclip/venv/environment.sh`.
  
- For all others:
    - Install Python (3.6+)
    - Install peak (`pip install eclip-peak`)
    - Install [IDR](https://github.com/nboley/idr) (2.0.3+)
    - Install Perl (5.10.1+) with the following packages:
        - Statistics::Basic (`cpanm install Statistics::Basic`)
        - Statistics::Distributions (`cpanm install Statistics::Distributions`)
        - install Statistics::R (`cpanm install Statistics::R`)
    
## Usage
- For Van Nostrand Lab
  
    After activate peak's environment call` peak -h` to see the detailed usage. 


- For all others:

    After successfully installed Python, peak, Perl (with required packages), 
    call `peak -h` inside your terminal to see the following detailed usage:
  
```shell
$ peak -h
usage: peak [-h] 
            [--ip_bams IP_BAMS [IP_BAMS ...]] 
            [--input_bams INPUT_BAMS [INPUT_BAMS ...]] 
            [--peak_beds PEAK_BEDS [PEAK_BEDS ...]] 
            [--read_type READ_TYPE] [--outdir OUTDIR] 
            [--species SPECIES] 
            [--l2fc L2FC] [--l10p L10P] [--idr IDR] 
            [--dry_run] [--cores] [--debug]

Pipeline for using IDR to identify a set of reproducible peaks given eClIP dataset 
with two or three replicates.

optional arguments:
  -h, --help            show this help message and exit
  --ip_bams IP_BAMS [IP_BAMS ...]
                        Space separated IP bam files (at least 2 files).
  --input_bams INPUT_BAMS [INPUT_BAMS ...]
                        Space separated INPUT bam files (at least 2 files).
  --peak_beds PEAK_BEDS [PEAK_BEDS ...]
                        Space separated peak bed files (at least 2 files).
  --ids IDS [IDS ...]   Optional space separated short IDs (e.g., S1, S2, S3) for datasets.
  --read_type READ_TYPE
                        Read type of eCLIP experiment, either SE or PE.
  --outdir OUTDIR       Path to output directory.
  --species SPECIES     Short code for species, e.g., hg19, mm10.
  --l2fc L2FC           Only consider peaks at or above this l2fc cutoff, default: 3.
  --l10p L10P           Only consider peaks at or above this l10p cutoff, default: 3.
  --idr IDR             Only consider peaks at or above this idr score cutoff, default: 0.01.
  --cores CORES         Maximum number of CPU cores for parallel processing, default: 1.
  --dry_run             Print out steps and inputs/outputs of each step without 
                        actually running the pipeline.
  --debug               Invoke debug mode (only for develop purpose).

```
  
## Outline of workflow
 - Normalize CLIP IP BAM over INPUT for each replicate
 - Peak compression/merging on input-normalized peaks for each replicate
 - Entropy calculation on IP and INPUT read probabilities within each peak for each replicate
 - Run IDR on peaks ranked by entropy
 - Normalize IP BAM over INPUT using new IDR peak regions
 - Identify reproducible peaks within IDR regions

## Examples

- eCLIP with 2 replicates
    
    Assuming we have eCLIP pipeline run successfully and have the following files generated 
    for species `hg19`:
    ```
    replicate 1:
        IP BAM: ip1.bam
        INPUT BAM: input1.bam
        Peak BED: clip1.peak.clusters.bed
    replicate 2:
        IP BAM: ip2.bam
        INPUT BAM: input2.bam
        Peak BED: clip2.peak.clusters.bed
    ```
  
    The pipeline then can be called like this to identify reproducible peaks:
    ```shell
    peak \
        --ip_bams ip1.bam ip2.bam \
        --input_bams input1.bam input2.bam \
        --peak_beds clip1.peak.clusters.bed clip2.peak.clusters.bed \
        --species hg19
    ```
  
- eCLIP with 3 replicates
    
    Assuming we have eCLIP pipeline run successfully and have the following files generated 
    for species `hg19`:
    ```
    replicate 1:
        IP BAM: ip1.bam
        INPUT BAM: input1.bam
        Peak BED: clip1.peak.clusters.bed
    replicate 2:
        IP BAM: ip2.bam
        INPUT BAM: input2.bam
        Peak BED: clip2.peak.clusters.bed
    replicate 3:
        IP BAM: ip3.bam
        INPUT BAM: input3.bam
        Peak BED: clip3.peak.clusters.bed
    ```
  
    The pipeline then can be called like this to identify reproducible peaks:
    ```shell
    peak \
        --ip_bams ip1.bam ip2.bam ip3.bam \
        --input_bams input1.bam input2.bam input3.bam \
        --peak_beds clip1.peak.clusters.bed clip2.peak.clusters.bed clip3.peak.clusters.bed \
        --species hg19
    ```
Note:

 - The indentation of the command does not matter, you can write it on the same line.
 - The order of bam and peak files followed by `--ip_bams`, `input_bams`, and `peak_beds` 
   DOES matter, make sure you pass them in a consistent order for these three parameters.
 - There are 3 cutoffs can be set for fine tune the peak filtering, see Usage part for 
   more details.
 - If the pipeline failed, check the log to identify the error and make necessary changes,
   re-run the pipeline will skip successfully processed parts only continue to processed 
   failed and unprocessed parts.
   
## Output
The peak pipeline will output 5 different types of files into the current work directory 
or into a user specified output directory (via `--outdir`):
1. *.bed: either a 6 columns or 9 columns bed file saves information for peaks.
2. *.tsv: TSV separated text file saves more information in addition to the BED file.
3. *.txt: text file saves the mapped reads count
4. *.out: TAB separated text file generated by IDR.
5. *.png: plot generated by IDR.

All filenames of output files are self-explained, only the basename of peak bed files (
after the removal of .peak.clusters.bed) was used to mark the name of each replicate.

The reproducible peaks can be found in 
*.reproducible.peaks.bed and additional information can be found in *.reproducible.peaks.custom.tsv.
While the former file is 6-column bed file, the later one is a TSV separated text file with the 
following columns in order:
- IDR region (entire IDR identified reproducible region)
- Peak (reproducible peak region)
- Geomean of the l2fc
- Columns of log2 fold change (2 or 3 columns for 2 or 3 replicates experiment, respectively)
- Columns of -log10 p-value (2 or 3 columns for 2 or 3 replicates experiment, respectively)
