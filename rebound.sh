
#!/bin/bash
#.
#.                         |__  o\
#.                         | W    \O
#.      ###############    |       |\_           |\
#.      rebound.sh v1.0    |      /-\            \O
#.      ###############    |    /     \           |
#.                         |                     /|
#.                         |                    |  \
#.
#. Usage:
#.        rebound.sh  [ -r ] [ FASTA ]
#.                    [ -o ]
#.                    [ -m ]
#.                    [ -s ]
#.                    [ -c ]
#.                    [ -d ]
#.                    [ -t ]  BAM1, BAM1, ...
#.
#. Description:
#.      This script adds SLAM-seq tags to BAM files, thus allowing to use
#.      any aligner of choice for downstream analysis with the SlamDunk
#.      package.
#.
#. Requirements:
#.      Samtools    ( https://anaconda.org/bioconda/samtools )
#.      BamUtils    ( https://anaconda.org/bioconda/bamutil )
#.      GAWK        ( https://anaconda.org/anaconda/gawk )
#.
#. Options:
#.  -r   fasta reference file; should be the reference file used
#.       to generate the BAM file(s) and have .fa or .fasta extension
#.  -o   output directory name (saves in working directory if not set)
#.  -m   add MD tags to BAM files using samtools calmd; will overwrite
#.       MD tag if such already exist, and add/modify NM tag accordingly;
#.       requires -r flag set
#.  -p   treat input BAM files as paired-end (default is single-end)
#.  -c   skip clipping of paired reads overlaps in paired-end data
#.       (not recommended); by default one read from every
#.       pair is clipped to avoid double counting of T->C conversions
#.       within the overlapping region; requires -r flag set
#.  -t   number of CPU threads to use (affects samtools)
#.  -d   whether to remove duplicate reads based on their bitwise flag
#.       (requires BAM files to have duplicates marked with Pickard
#.       MarkDuplicates or similar)
#.  BAM  (required) list of BAM files to process (wildcards * accepted)
#.  -h   prints this help

help() {
grep "#\." $0 | sed 's/#\.//'
exit 0
}

usage() {
   echo "Usage: rebound.sh   [ -r ] [ FASTA ]
                    [ -o ]
                    [ -m ]
                    [ -p ]
                    [ -c ]
                    [ -d ]
                    [ -t ]  BAM1, BAM1, ...

To get more help type: rebound.sh -h"
}

exit_err() {
  usage
  exit 1
}

[ $# -eq 0 ] && exit_err

# Set defaults
THR_ARG=1
OUTDIR_ARG="."
while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts ":o:r:mpct:dh" opcje
    do
    case $opcje in
        o) # Specify a desctination folder for the results.
            OUTDIR_ARG=${OPTARG}
            ;;
        r) # Specify a valid reference file.
            REF_ARG=${OPTARG}
            ;;
        m) # Whether to overwrite MD tags.
            MD_ARG=1
            ;;
        p) # Toogle sinle-end or paired-end (default is single-end).
            PAIRED_ARG=1
            ;;
        c) # Wether to filter out dups.
            NODUPS_ARG=1
            ;;
        c) # Wether to clip paird reads overlaps.
            NOCLIP_ARG=1
            ;;
        t) # Number of threads for samtools.
            THR_ARG=${OPTARG}
            ;;
        h) # Display help.
            help
            ;;
        esac
    done
    shift "$((OPTIND-1))"
    ARGS="${ARGS} $1"
    shift
done

BAM_LIST=$ARGS

check_soft() {
    if ! command -v samtools &> /dev/null ; then
       echo "Error. Samtools not found."
       exit 1
   elif ! command -v bam &> /dev/null ; then
       echo "Error. BamUtils not found."
       exit 1
   else
       echo Found Samtools $(samtools --version | head -1 | cut -d' ' -f2)
       echo Found BamUtils $(bam clipOverlap 2>&1 | head -1 | cut -d' ' -f2 | sed 's/;//g')
   fi
}

check_soft

##### parameter check #####

if [[ $MD_ARG == 1 ]]; then
    if ! [[ $REF_ARG =~ ^.+\.fasta$  ||  $REF_ARG =~ ^.+\.fa$ ]]; then
        echo "[ -m ] flag requires a valid reference file and $REF_ARG does not look like one..."
        exit_err
    else
        echo "MD tags will be added to reads and pre-existing MD tags will be overwritten"
    fi
fi

if ! [[ $THR_ARG =~ ^[0-9]+$ ]]; then
    echo "Threads # must be an integer!"
  exit_err
fi

if [[ $PAIRED_ARG == 1 ]]; then
   echo "Treat input as paired-end is ON"
else
   echo "Treat input as paired-end is OFF"
fi

if [[ $PAIRED_ARG == 1 ]]; then
    if [[ $PAIRED_ARG == 1 ]]; then
       echo "Treat input as paired-end is ON"
    else
   echo "Clipping of paired read overlaps will NOT be clipped"
else
   if ! [[ $REF_ARG =~ ^.+\.fasta$  ||  $REF_ARG =~ ^.+\.fa$ ]]; then
       echo "[ -c ] flag requires a valid reference file and $REF_ARG does not look like one"
       exit_err
   fi
   echo "Soft-clipping of the overlapping portions of read pairs is ON"
fi

if [[ $BAM_LIST == " " ]]; then
     echo "No BAM file(s) provided!"
   exit_err
else
    for BAM in $BAM_LIST
    do
        if ! [[ $BAM =~ ^.+\.bam$ ]] ; then
            echo $BAM "is not a valid file name!"
            exit_err
        fi
    done
fi

echo "$REF_ARG is used as a reference file"
echo "Using $THR_ARG CPU thread(s)"

if [[ $OUTDIR_ARG == "." ]]; then
    echo "Output files will be saved in the working directory as the [ -o ] flag was not set"
else
    mkdir -p $OUTDIR_ARG
    mkdir -p $OUTDIR_ARG/stderr
fi

##### main routine #####

for BAM in $BAM_LIST
do
    BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' )
    BAM_FILE=$( echo $BAM | xargs basename )

    echo "Processing $BAM_FILE"
    echo "Saving Samtools stderr output to $OUTDIR_ARG/stderr/$BAM_NAME\_st_stderr.txt"

    if [[ $NOCLIP_ARG != 1 ]]; then
        echo Saving ClipOverlap stderr output to $OUTDIR_ARG/stderr/$BAM_NAME\_clip_stderr.txt
    fi
    echo "Preparing SLAMseq tags"
    if [[ $NODUPS_ARG == 1 ]]; then
        samtools view -@ $THR_ARG -h -F 1024 $BAM 2> $OUTDIR_ARG/stderr/$BAM_NAME\_st_stderr.txt
    else
        samtools view -@ $THR_ARG -h $BAM 2> $OUTDIR_ARG/stderr/$BAM_NAME\_st_stderr.txt
    fi |\

    if [[ $MD_ARG == 1 ]]; then
        samtools calmd -@ $THR_ARG - $REF_ARG 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_st_stderr.txt
    fi |\

    if [[ $NOCLIP_ARG != 1 ]]; then
        samtools sort -n -@ $THR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_st_stderr.txt |\
        bam clipOverlap --in - --out - --storeOrig CG --readName --stats 2> $OUTDIR_ARG/stderr/$BAM_NAME\_clip_stderr.txt
    fi |\
    gawk -f rebound.awk 2> $OUTDIR_ARG/stderr/$BAM_NAME\_gawk_stderr.txt 1> $OUTDIR_ARG/$BAM_NAME\_rebound.sam
echo "GAWK stderr output saved to $OUTDIR_ARG/stderr/$BAM_NAME\_gawk_stderr.txt"
echo "Output BAM saved to $OUTDIR_ARG/$BAM_NAME\_rebound.sam"
done
printf "\nAll done!\n"
