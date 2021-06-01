                         |__  o\
                         | W    \O
      ###############    |       |\_           |\
      rebound.sh v1.0    |      /-\            \O
      ###############    |    /     \           |
                         |                     /|
                         |                    |  \



**Intended use:**
* Quantifying T->C conversions in RNA sequencing data
* Processing paired-end BAM files for downstream analysis (quantification of T->C conversions) with SlamDunk package (https://anaconda.org/bioconda/slamdunk)

**How to run:**
* Make script executable (`chmod u+x rebound.sh`) and run (see examples). Because GAWK is single-threaded, it is best to run the script using only one thread (flag `-t 1`). For multiple file script can be run in parallel, as shown in one of the examples below.

**Input:**
Paired-end or single-end SAM/BAM files, annotation file (gtf/gff) and genomic reference.

**Output:**  
* Table of counts of mismatches across all possible nucleotide combinations
* Table of counts of mismatches across all possible nucleotide combinations (background-corrected)
* Sanity-check plots

**Prerequisites:**
* Samtools    (https://anaconda.org/bioconda/samtools)
* BamUtils    (https://anaconda.org/bioconda/bamutil)
* GAWK        (https://anaconda.org/anaconda/gawk)<br>

**Options and usage:**<br>

`rebound.sh [options] [BAMs]`

| Option      | Default     | Description |
| ----------- | ----------- | ----------- |
|`-a`|`./`|Path to `rebound.awk`|
|`-o`|`./`|Output directory|
|`-l`|`star`|Aligner used to prepare BAM files (one of `bbmap` or `star`). Rebound was tested with output of STAR and BBmap (sam=1.3 flag). If using other aligner, alignment files should be in SAM 1.3 (CIGARs must not contain X and = operations, but only MNDIS). Files in SAM 1.4 format can be converted to SAM 1.3 with `reformat.sh in=input.bam out=output.bam sam=1.3`.|
|`-r`||Fasta reference file with full path. Genome reference that is the same as the one used to prepare BAM files or updated based on SNP analysis (variant calling). The latter can be done with `bcftools` (see examples).|
|`-f`||Overwrite MD tags. With this flag set MD tags are added if absent, or overwritten if present. In mode "normal" (default) MD tag are replaced by default, so this flag has no effect.|
|`-i`|`0`|Filter reads containing indels smaller than this value (reads will be excluded from the output).|
|`-m`|`normal`| Mode of function (`slamdunk` or `normal`). In mode `normal` a table with mismatches per gene and relevant sanity check plots are produced. In mode `slamdunk` a slamDunk-compatible BAM file(s) are produced in addition. Currently, mode `slamdunk` works best with single-end input but if paired-end BAM files are provided (type of file should be determined automatically), rebound will attempt to convert them to single-end while retaining original mapping positions (experimental).|
|`-g`||Annotation (gtf/gff) file with path. If using unstranded libraries, this file may include strandedness information appended to gene identifiers as double underscore followed by `+` or `-` sign, e.g., `ENSMUSG00000079037_+`.|
|`-q`|`254,20,0.95,27,27` (star) <br>`1,100000,0.95,27,27` (bbmap)| Quality filtering settings. It is a string of 5 comma-separated values (order of values important) defining: 1. MAPQ (integer): mapping quality filter; the default value for star (255) keeps only reads with `MAPQ 255`, which is STAR's default value for unique mappers. 2. ED (integer): Read filtering threshold by edit distance (e.g., 20 keeps reads with values below this threshold). Set to large values to disable the filter (this may be necessary if introns (N) are represented as deletions (D) in CIGARs) 3. XF (decimal [0-1): Read filtering threshold by mismatch fraction (e.g., 0.95 excludes reads with mismatch content above 5%). 4. QV (integer [1-41]): Base quality threshold; excludes mismatches when called base quality is below this value (range [1-41]). 5. QVTC (range [1-41]): Same as QV but only affects T->C mismatches. Must be higher then QV to have an effect.|
|`-s`|`forward`|Strandedness. In case of unstranded libraries Rebound can read strand information from annotation file. In such case strandedness info (`+` or `-`) must be appended to gene identifiers in the GTF file (double underscore followed by `+` or `-` sign, e.g., `ENSMUSG00000079037__+`). GTF file can be modified as shown in the examples.|
|`-d`||Remove duplicated reads. Requires duplicates marked with Picard MarkDuplicates or a similar tool. If this flag is NOT set, and duplicates are marked, this information will be included in the final output counts, allowing two distinct analyses (with and without duplicate removal).|
|`-u`||Use only uniqly mapped reads. Default is to use all reads.|
|`-x`|`1G`|Per thread memory block size (affects samtools).|
|`-t`|`1`|Number of threads (affects samtools).|
|`-h`||Prints help, overrides any remaining options.|
|`BAM`|List BAMs to process (wildcards [*] accepted). Also accepts SAM format.|

**Examples**    

Adding strand information to GTF file (works on GENCODE GFT format):

`cat annotations.gtf |
awk 'BEGIN{OFS="\t"}{gsub(/\.[[:digit:]]+";|";/,"__"$7"\";",$10); print $0}' > annotations_with_strand.gtf`

Running rebound on a STAR-aligned BAM file, reverse strandedness, keep unique reads only, ignore reads with deletions <= 30:  

`rebound.sh -r genome.fa -o output/dir -l star -g annotations.gtf -s reverse -u -f -i 30 my.bam`

**Output description**  
![Summary table](/img/read_stats.jpg)
Counts of reads passing filters and a summary of filter settings.
* *Total reads*: all reads initially considered.
* *Reads passing indel filter*: reads passing indel filter.
* *Reads with matches/mismatches* - all reads reads with M operations in CIGARs (reads with no M operations may occur due to clipping of overlapping pairs in paired-end data).
* *Reads with mismatches*: fraction that contain mismatched bases.
* *Reads with mismatches filtered* - fraction of reads with mismatches that passed QC filters [ -q ].
* *Reads with mismatches filtered (unique assignment)*:  (mode normal only) fraction of reads with mismatches, passing QC filters, and unambiguously assigned to a feature in GTF file. Specifically, it excludes all entries marked by HTseq as any of the following: `__ambiguous|__no_feature|_too_low_aQual|__not_aligned|__alignment_not_unique`. The count includes multi-mappers if  [ -u ] flag was not set.
* *Reads without mismatches*: fraction of reads that contain mismatched bases.
* *Reads without mismatches filtered*: fraction of reads without mismatches that passed QC filters [ -q ].
* *Reads without mismatches filtered (unique assignment)*: (mode normal only) fraction of reads without mismatches, passing QC filters, and unambiguously assigned to a feature in GTF file. Specifically, it excludes all entries marked by HTseq as any of the following: `__ambiguous|__no_feature|_too_low_aQual|__not_aligned|__alignment_not_unique`. The count includes multi-mappers if [ -u ] flag was not set.

![Rates main](/img/coverage_histogram.png)
Histogram of coverage over tymines (T). Red, blue and black lines denote 100, 1000, and 10000 cutoffs respectively. Genes and ERCC spike-ins are automatically discriminated provided that Ensembl gene identifiers are used, and ERCCs are included in the GTF file.

![Rates plot](/img/rates_genes_plot.png)
Mismatch rates calculated (**for each gene independently**) by dividing the number of observed mismatches by the total count of reference nucleotides of the corresponding type (e.g., T for T->C mismatch). <br>
![Rate equation](/img/rate_equation.png)

![Rates coverage plot](/img/rates_per_cov_genes_plot.png)
Same as the previous plot, but including coverage bins.

![Overrepresentation plot](/img/overrepresentation_plot_bulk.png)
Overrepresentation or T->C mismatches across the landscape of all possible (ACTG) mismatches (left facet) or only T mismatches (right facet). Useful to determine optimal read filter parameters to improve signal-to-noise ratio for SLAM-seq analysis.

![Overrepresentation plot with coverage](/img/overrepresentation_plot_w_coverage.png)
Same as the previous plot, but including coverage bins.

![Overrepresentation plot with coverage](/img/coverage_overrepresentation_histogram.png)
Same data as the previous plot, but shown as histograms.

![Overrepresentation plot with coverage](/img/rates_spike_plot.png)
Mismatch rates for the optional in-vitro transcribed control transcript (slam). This transcript serves as internal control of alkylation efficacy and can be used to normalize for potential variation between samples. <br> "Slam" is an optional in-vitro transcribed mRNA containing synthetic 4-thiouridines only (no "normal" uridines). It serves as a control for alkylation efficacy and may be used to normalize for variation between samples. "Slam" spike-in control is typically added to RNA after initial purification and prior to alkylation.

**References**<br>
Herzog, V., Reichholf, B., Neumann, T. et al. Thiol-linked alkylation of RNA to assess expression dynamics. *Nat Methods* 14, 1198–1204 (2017). https://doi.org/10.1038/nmeth.4435
