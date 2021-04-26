
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
#.        rebound.sh  [ -a ]
#.                [ -o ]
#.                [ -r ] [ FASTA ]
#.                [ -m ]
#.                [ -p ]
#.                [ -c ]
#.                [ -n ] (do not use)
#.                [ -x ] (do not use)
#.                [ -e ] (do not use)
#.                [ -s ] (must be rev)
#.                [ -d ]
#.                [ -d ]
#.                [ -t ]  BAM1, BAM1, ...
#.
#. Description:
#.      This script adds SLAM-seq tags to BAM files, and thus allows to use
#.      STAR or other aligner of choice for downstream analysis with the slamDunk
#.      package. Can be used downstream of nf-core/rna-seq pipeline.
#.      Currently works only for stranded sequencing with "reverse" strandedness.
#.
#. Prerequisites:
#.      Samtools    ( https://anaconda.org/bioconda/samtools )
#.      BamUtils    ( https://anaconda.org/bioconda/bamutil )
#.      GAWK        ( https://anaconda.org/anaconda/gawk )
#.      slamDunk    ( https://anaconda.org/bioconda/slamdunk )
#.
#. Options:
#. -a (./)  Path to rebound.awk script.
#. -o (./)  Output directory.
#. -r       Fasta reference file. Should be the same file that was used
#.          to generate the BAM file(s) and have .fa or .fasta extension
#. -m       Force replacement of N with D in CIGARs and generations of MD tags
#.          using samtools calmd. Will overwrite MD tag if already present, and
#.          add/modify NM tag accordingly. Requires -r flag set. This flag is set
#.          automatically if N operations are present in any of the CIGARs
#.          (to test for that the first 10000 BAM records are probed).
#. -p       Treat input BAM files as paired-end (default is single-end).
#. -c       Clip paired-end reads overlaps in paired-end data (not recommended).
#.          By default one read from every pair is clipped to avoid double
#.          counting of T->C conversions within the overlapping region.
#.          Requires -r flag set. After clipping, samtools calmd is run
#.          to update MD and NM tags.
#. -n (254) Read filtering threshold by MAPQ. Default (254) keeps only reads
#.          with MAPQ 255, which is STAR's default value for unique mappers.
#.          Uses slamDunk filter. (to be implemented)
#. -x (.95) Read filtering threshold by mismatch fraction (.95 preserves
#.          reads with at least 95% of matching bases). Uses slamDunk filter. (to be implemented)
#. -t (1)   Number (integer) of CPU threads to use. Affects samtools.
#. -d       Force removing of duplicate reads based on bitwise flags.
#.          Requires BAM files to have duplicates marked with Picard
#.          MarkDuplicates or a similar tool.
#. -g (1G)  Per thread memory block size for samtools sort.
#. BAM      (required) BAM files to process (lists and wildcards [*] accepted).
#. -h       Prints this help and overrides any remaining flags.

help() {
grep "#\." $0 | sed 's/#\.//'
exit 0
}

REBOUND_COMMAND="rebound.sh $@"

usage() {
   echo "Usage: rebound.sh   [ -a ]
                    [ -o ]
                    [ -r ] [ FASTA ]
                    [ -m ]
                    [ -p ]
                    [ -c ]
                    [ -n ] (do not use)
                    [ -x ] (do not use)
                    [ -e ] (do not use)
                    [ -s ] (must be rev)
                    [ -d ]
                    [ -g ]
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
OUTDIR_ARG="./"
AWK_FILE="./rebound.awk"

while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts ":a:o:r:mpcs:dg:t:h" opcje
    do
    case $opcje in
        a) # Specify a path to AWK script.
            AWK_FILE="$(echo ${OPTARG} | sed "s/\\/$//")/rebound.awk"
            ;;
        o) # Specify a destination folder for the results.
            OUTDIR_ARG=${OPTARG}
            ;;
        r) # Specify a valid reference file with full path.
            REF_ARG=${OPTARG}
            ;;
        m) # Whether to overwrite MD tags.
            MD_ARG=1
            ;;
        p) # Toogle sinle-end or paired-end (default is single-end).
            PAIRED_ARG=1
            ;;
        c) # Whether to clip paird reads overlaps.
            CLIP_ARG=1
            ;;
        s) # Strandedness.
            STRANDEDNESS_ARG=${OPTARG}
            ;;
        d) # Whether to filter out dups.
            NODUPS_ARG=1
            ;;
        g) # Memory for samtools sort in G.
            MEM_ARG=${OPTARG}
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
    if ! [ -f "$AWK_FILE" ]; then
        echo "Fatal error: $AWK_FILE not found."
        exit 1
    elif ! command -v samtools &> /dev/null ; then
        echo "Fatal error. Samtools not found."
        exit 1
    elif ! command -v bam &> /dev/null ; then
        echo "Fatal error: BamUtils not found."
        exit 1
    else
        echo Found Samtools $(samtools --version | head -1 | cut -d' ' -f2)
        echo Found BamUtils $(bam clipOverlap 2>&1 | head -1 | cut -d' ' -f2 | sed 's/;//g')
   fi
}

check_soft

##### parameter check #####

if [[ $MD_ARG == 1 ]]; then
    if ! [ -f "$REF_ARG" ] ; then
        echo "Fatal error: $REF_ARG not found."
        exit 1
    elif ! [[ $REF_ARG =~ ^.+\.fasta$  ||  $REF_ARG =~ ^.+\.fa$ ]] ; then
        echo "[ -m ] flag requires a valid reference file and $REF_ARG does not look like one..."
        exit_err
    else
        echo "Fill MD tags [ -m ]: ON"
    fi
fi

if ! [[ $THR_ARG =~ ^[0-9]+$ ]] ; then
    echo "Threads # must be an integer!"
    exit_err
fi

if [[ $CLIP_ARG == 1 ]]; then
   if [[ $PAIRED_ARG == 1 ]]; then
       if ! [[ $REF_ARG =~ ^.+\.fasta$ || $REF_ARG =~ ^.+\.fa$ ]] ; then
           echo "Error: [ -c ] flag requires a valid reference file and $REF_ARG does not look like one"
           exit_err
       fi
       echo "Paired-end mode [ -p ]: ON"
       echo "Clip overlaps [ -c ]: ON"
   else
       echo "Paired-end mode [ -p ]: OFF"
       CLIP_ARG=0
       echo "Clip overlaps [ -c ] automatically re-set to OFF"
   fi
else
    if [[ $PAIRED_ARG == 1 ]]; then
        echo "Paired-end mode [ -p ]: ON"
        echo "Clip overlaps [ -c ]: OFF"
    else
        echo "Paired-end mode [ -p ]: OFF"
    fi
fi

if [[ $BAM_LIST == " " ]]; then
     echo "No BAM file(s) provided!"
   exit_err
else
    for BAM in $BAM_LIST
    do
        if ! [[ $BAM =~ ^.+\.bam$ ]] ; then
            echo $BAM "is not a valid input file name!"
            exit_err
        fi
    done
fi

echo "$REF_ARG is used as genome reference"
echo "Using $THR_ARG CPU thread(s)"
echo "Using $MEM_ARG RAM per thread"

if [[ $OUTDIR_ARG == "." ]]; then
    echo "Saving output files in the current directory as the [ -o ] flag was not set"
else
    mkdir -p $OUTDIR_ARG
    mkdir -p $OUTDIR_ARG/stderr
fi

##### main routine #####
START_TIME=$SECONDS
for BAM in $BAM_LIST
do
    BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' )
    BAM_FILE=$( echo $BAM | xargs basename)

    printf " ------------------\n Processing %s\n ------------------\n" "$BAM_FILE"

    samtools view -@ $THR_ARG $BAM | head -10000 | gawk '{if ($6 ~ /N/){exit 0}}'

    if [ $? -eq 0 ]; then
        echo "Found Ns in CIGARs, will be replaced with Ds"
        echo "(MD tags will be automatically updated)"
        NINCIG=1
    fi

    echo "Strandedness is set to $STRANDEDNESS_ARG"
    echo "Preparing SLAM-seq tags"

    if [[ $NODUPS_ARG == 1 ]]; then
        samtools view -@ $THR_ARG -h -F 1024 $BAM 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_samtools_stderr.txt
    else
        samtools view -@ $THR_ARG -h $BAM 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_samtools_stderr.txt
    fi |\

    if [[ $MD_ARG == 1 ]]; then
        if [[ $(samtools view -H $BAM | head -1) =~ ^.*oordinate.*$ ]]; then
            cat
        else
            samtools sort -m $MEM_ARG -@ $THR_ARG -T $OUTDIR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_samtools_stderr.txt
        fi |\
        samtools calmd -@ $THR_ARG - $REF_ARG 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_calmd_stderr.txt
    else
        cat
    fi |\

    if [[ $CLIP_ARG == 1 ]]; then
        samtools sort -n -m $MEM_ARG -@ $THR_ARG -T $OUTDIR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_samtools_stderr.txt |\
        bam clipOverlap --in - --out - --readName --storeOrig OC --stats 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_clipOverlap_stderr.txt |\
        samtools sort -m $MEM_ARG -@ $THR_ARG -T $OUTDIR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_samtools_stderr.txt
    else
        cat
    fi |\

    if [[ $NINCIG == 1 ]]; then
        gawk 'BEGIN{OFS="\t"}; function tag_finder(tag, i){for (i=1; i<=NF+1; i++){if (i==NF+1) {print "Fatal error: unable to find " tag " at record " $0; exit 1 } else if (substr($i,1,5)!=tag) {continue} else { return $i }}} {if ($1 ~ /^@/) { print $0 } else { gsub($6, gensub("N", "D", "g", $6), $0); if ($0 ~ /^.*NM:i:.*$/) { print $0, "ED:i:" gensub("NM:i:", "", "g", tag_finder("NM:i:")) } else { print $0 }}}' |\
        samtools calmd -@ $THR_ARG - $REF_ARG 2>> $OUTDIR_ARG/stderr/$BAM_NAME\_calmd_stderr.txt
    else
        cat
    fi |\

    gawk -v REBOUND_COMMAND="$REBOUND_COMMAND" -v BAM_NAME="$BAM_NAME" -v STRANDEDNESS="$STRANDEDNESS_ARG" -f $AWK_FILE 2> $OUTDIR_ARG/stderr/$BAM_NAME\_rebound_stderr.txt 1> $OUTDIR_ARG/$BAM_NAME\_reboundbig.sam

    if [ $? -eq 0 ]; then
        echo "Processing $BAM_FILE completed successfully"
        echo "Output saved to $OUTDIR_ARG/${BAM_NAME}_rebound.sam"
    else
        echo "Processing $BAM_FILE FAILED"
    fi

    if [ $(stat $OUTDIR_ARG/stderr/$BAM_NAME\_calmd_stderr.txt | awk '{print $1}') -gt 100000 ]; then
        tail -1000 $OUTDIR_ARG/stderr/$BAM_NAME\_calmd_stderr.txt > $OUTDIR_ARG/stderr/calmd_stderr.tmp && cat $OUTDIR_ARG/stderr/calmd_stderr.tmp > $OUTDIR_ARG/stderr/$BAM_NAME\_calmd_stderr.txt
    fi

    [ -f $OUTDIR_ARG/tmp_rebound.bam ] && rm $OUTDIR_ARG/tmp_rebound.bam
    [ -f $OUTDIR_ARG/stderr/calmd_stderr.tmp ] && rm $OUTDIR_ARG/stderr/calmd_stderr.tmp
done

TOTAL_TIME=$(($SECONDS - $START_TIME))
echo "Finished processing" $(echo "$BAM_LIST" | wc -w | awk '{print $1}') "file(s) in" $(echo "scale=2;$TOTAL_TIME/3600" | bc) "hour(s)."
