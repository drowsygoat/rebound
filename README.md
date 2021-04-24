                         |__  o\
                         | W    \O
                         |       |\_           |\
      rebound.sh v1.0    |      /-\            \O
                         |    /     \           |
                         |                     /|
                         |                    |  \

     Usage:
            rebound.sh  [ -a ]
                        [ -o ]
                        [ -r ] [ FASTA ]
                        [ -m ]
                        [ -p ]
                        [ -c ]
                        [ -n ]
                        [ -x ]
                        [ -t ]
                        [ -d ]  BAM1, BAM1, ...

     Description:
          This script adds SLAM-seq tags to BAM files, and thus allows to use
          STAR or other aligner of choice for downstream analysis with the slamDunk
          package. Can be used downstream of nf-core/rna-seq pipeline.

     Prerequisites:
          Samtools    ( https://anaconda.org/bioconda/samtools )
          BamUtils    ( https://anaconda.org/bioconda/bamutil )
          GAWK        ( https://anaconda.org/anaconda/gawk )
          slamDunk    ( https://anaconda.org/bioconda/slamdunk )

     Options:
     -a (./)  Path to rebound.awk script.
     -o (./)  Output directory.
     -r       Fasta reference file. Should be the same file that was used
              to generate the BAM file(s) and have .fa or .fasta extension
     -m       Force replacement of N with D in CIGARs and generations of MD tags
              using samtools calmd. Will overwrite MD tag if already present, and
              add/modify NM tag accordingly. Requires -r flag set. This flag is set
              automatically if N operations are present in any of the CIGARs
              (to test for that the first 10000 BAM records are probed).
     -p       Treat input BAM files as paired-end (default is single-end).
     -c       Clip paired-end reads overlaps in paired-end data (not recommended).
              By default one read from every pair is clipped to avoid double
              counting of T->C conversions within the overlapping region.
              Requires -r flag set. After clipping, samtools calmd is run
              to update MD and NM tags.
     -n (254) Read filtering threshold by MAPQ. Default (254) keeps only reads
              with MAPQ 255, which is STAR's default value for unique mappers.
              Uses slamDunk filter. (to be implemented)
     -x (.95) Read filtering threshold by mismatch fraction (.95 preserves
              reads with at least 95% of matching bases). Uses slamDunk filter. (to be implemented)
     -t (1)   Number (integer) of CPU threads to use. Affects samtools.
     -d       Force removing of duplicate reads based on bitwise flags.
              Requires BAM files to have duplicates marked with Picard
              MarkDuplicates or a similar tool.
     BAM      (required) BAM files to process (lists and wildcards [*] accepted).
     -h       Prints this help and overrides any remaining flags.
