
#!/bin/bash
#.
#.                         |__  o\
#.                         | W    \O
#.      ###############    |       |\_           |\
#.      rebound.sh v1.1    |      /-\            \O
#.      ###############    |    /     \           |
#.                         |                     /|
#.                         |                    |  \
#.
#. Usage:
#.        rebound.sh [OPTIONS] [BAM FILES]
#.
#. Prerequisites:
#.      Samtools    ( https://anaconda.org/bioconda/samtools )
#.      BamUtils    ( https://anaconda.org/bioconda/bamutil )
#.      GAWK        ( https://anaconda.org/anaconda/gawk )
#.
#. Options:
#. -a (./)  Path to rebound.awk script.
#. -o (./)  Output directory.
#. -m       (normal) Mode of function, one of <normal/slamdunk>
#. -r       Fasta reference file. This should be the file that was used to generate
#.          the BAM file(s) and have .fa or .fasta extension.
#. -m       Force replacement of N with D in CIGARs and generations of MD tags
#.          using samtools calmd. Will overwrite MD tag if already present, and
#.          add/modify NM tag accordingly. Requires -r flag set. This flag is set
#.          automatically if N operations are present in any of the CIGARs
#.          (to test for that the first 10000 BAM records are probed).
#. -p       <pe/se> Whether the input is paired-end (pe) or single-end (se) data (se).
#. -s       <forward/reverse> Strandedness of sequencing data (SLAM-seq requires
#.          stranded reads).
#. -q       <MAPQ(255),ED(10),MF(0.95),PHRED(27)>
#.          Read filtering settings:
#.          1. MAPQ/mapping quelity filter. Default (255) keeps only reads with MAPQ 255,
#.          which is STAR's default
#.          value for unique mappers.
#.          2. Read filtering threshold by edit distance (preserved reads with lower values).
#.          3. Read filtering threshold by mismatch fraction (.95 preserves reads with
#.          at least 95% of MATCHING bases).
#.          4. Ignore mismatches with called base quelity below this number (1-41.)
#.
#. -t       (1) Integer specifying the number of CPU threads to use. This setting only
#.          affects samtools. On a computing cluster, is is reccommended to run
#.          one instance of Rebound per sample using one thread for each instance.
#. -d       Force removal of duplicate reads based on bitwise flag.
#.          Requires BAM files to have duplicates marked with Picard
#.          MarkDuplicates or a similar tool. This is redundant in "normal" mode
#.          with flag [ -u ] set, as HTseq will not count such reads.
#. -u       Whether to count conversions using only uniqly mapped reads.
#. -g       (1G) Per thread memory block size for samtools sort.
#. BAM      BAM files to process (lists and wildcards [*] accepted).
#. -h       Prints this help and overrides any remaining flags.
# Options:
#   -o <path>
#   -a <file rebound.awk>
#   -r <file.fasta>
#   -g <file.gtf>
#   -q <INT,INT,DECIMAL,INT>
#   -s <reverse|forward>
#   -d
#   -u
#   -t INT
#   -x STRING
#   -h
#   "

NORM=$(tput sgr0)
BOLD=$(tput bold)
REV=$(tput smso)

help() {
grep "#\." $0 | sed 's/#\.//'
exit 0
}
REBOUND_COMMAND="rebound.sh $@"

usage() {
   echo "Usage:
   ${BOLD}rebound.sh [options] [file1.bam file2.bam ...]${NORM}
   Type ${BOLD}rebound.sh -h${NORM} to get more help."
}
exit_err() {
  usage
  exit 1
}
[ $# -eq 0 ] && exit_err
# Set some defaults
THR_ARG=1
MD_ARG=0
NODUPS_ARG=0
MODE_ARG="normal"
STRANDEDNESS_ARG="forward"
OUTDIR_ARG="./"
NINCIG=0
AWK_FILE="./rebound.awk"
while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts ":a:o:l:r:fi:g:m:q:s:dux:t:h" opcje
    do
    case $opcje in
        a) # Specify a path to AWK script.
            AWK_FILE="$(echo ${OPTARG} | sed "s/\\/$//")/rebound.awk"
            ;;
        o) # Specify a destination folder for the results.
            OUTDIR_ARG=${OPTARG}
            ;;
        l) # Specify a destination folder for the results.
            MAPPER_ARG=${OPTARG}
            ;;
        r) # Specify a valid reference file with full path.
            REF_ARG=${OPTARG}
            ;;
        f) # Force MD.
            MD_ARG=1
            ;;
        i) # Ban reads with indels.
            INDEL_BAN_ARG=${OPTARG}
            ;;
        g) # Specify a valid annotation (gtf) file with  path.
            GTF_ARG=${OPTARG}
            ;;
        m) # Mode of function.
            MODE_ARG=${OPTARG}
            ;;
        q) # Quality filtering settings (comma-separates numbers).
            QUALITY_ARG=${OPTARG}
            ;;
        s) # Strandedness.
            STRANDEDNESS_ARG=${OPTARG}
            ;;
        d) # Whether to filter out dups.
            NODUPS_ARG=1
            ;;
        u) # Include only unique reads.
            UNIQ_ARG=1
            ;;
        x) # Memory for samtools sort in G.
            MEM_ARG=${OPTARG}
            ;;
        t) # Number of threads for samtools.
            THR_ARG=${OPTARG}
            ;;
        h) # Display help.
            help
            ;;
        \?) #unrecognized option - show help
            echo -e \\n"Flag -${BOLD}$OPTARG${NORM} not allowed."
            usage
    esac
    done
    shift "$((OPTIND-1))"
    ARGS="${ARGS} $1"
    shift
done

BAM_LIST=$ARGS

##### check for installed software #####

check_soft() {
    if ! [ -f "$AWK_FILE" ]; then
        echo "Fatal error: $AWK_FILE not found."
        exit 1
    elif ! command -v samtools &> /dev/null ; then
        echo "Fatal error. Samtools could not be found."
        exit 1
    elif ! command -v bam &> /dev/null ; then
        echo "Fatal error: BamUtils could not be found."
        exit 1
    else
        echo Found Samtools $(samtools --version | head -1 | cut -d' ' -f2)
        echo Found BamUtils $(bam clipOverlap 2>&1 | head -1 | cut -d' ' -f2 | sed 's/;//g')
        if [[ $MODE_ARG == "slamdunk" ]]; then
            if ! command -v slamdunk -version &> /dev/null ; then
                echo "Fatal error: slamDunk could not be found."
                exit 1
            else
                echo Found $(slamdunk --v 2>&1 | head -1)
            fi
        fi
        if [[ $MODE_ARG == "normal" ]]; then
            if ! command -v htseq-count --help &> /dev/null ; then
                echo "Fatal error: HTseq could not be found."
                exit 1
            else
                echo "Found HTseq version" $(htseq-count --help | grep -o " \d\..*")
            fi
        fi
    fi
}

check_soft

##### parameter check #####
if ! [[ $MODE_ARG =~ ^.*slamdunk.*$|^.*normal.*$|^.*htseq_only.*$  ]]; then
    echo "Fatal error: Mode must be either \"slamdunk\" or \"normal\"" # or \"htseq_only\""
    exit_err
elif [[ $MODE_ARG =~ ^.*slamdunk.*$  ]] ; then
    echo "Mode: slamdunk"
elif [[ $MODE_ARG =~ ^.*htseq_only.*$  ]] ; then
    echo "Mode: htseq_only"
else
    echo "Mode: normal"
    if [[ $UNIQ_ARG == 1 ]]; then
       echo "Counting only unique reads [-u]"
    fi
fi

if [[ $OUTDIR_ARG == "." ]]; then
    echo "Output will be saved in the current directory"
else
    mkdir -p $OUTDIR_ARG
    mkdir -p $OUTDIR_ARG/stderr
fi

if [[ $MAPPER_ARG =~ ^.*bbmap.*$ ]] ; then
    echo "Aligner: BBmap"
elif [[ $MAPPER_ARG =~ ^.*star.*$ ]] ; then
    echo "Aligner: STAR"
fi

if [[ $MD_ARG -eq 1 ]] ; then
    echo "Forcing MD tags"
fi

samtools view -@ $((THR_ARG-1)) $(echo $BAM_LIST | cut -d" " -f1) | head -10000 | gawk '{if ($6 ~ /N/){exit 1}}'

if [[ $? -eq 1 ]]; then
    echo "${BOLD}Warning${NORM}: i(N)tron operations in CIGARs will be replaced with (D)eletion operations in the final output"
    NINCIG=1
else
    if [[ $INDEL_BAN_ARG =~ ban ]] ; then
        echo "${BOLD}Warning${NORM}: provided BAM file(s) seem to not discriminate introns from deletions and [ -i ] flag was set to \"ban\". This combination will effectively remove all reads that span splice junctions. In case this is not the expected behavior, set [ -i ] to a value in range of [1:30] instead."
    fi
fi

# indel filter
if ! [[ -z $INDEL_BAN_ARG ]] ; then
    if ! [[ ( $INDEL_BAN_ARG -gt 0 && $INDEL_BAN_ARG -lt 31 ) || $INDEL_BAN_ARG =~ ^.*disabled|ban.*$ ]] ; then
        echo "Indel threshold must be in range [1-30], \"disabled\" (default) or \"ban\"."
        exit_err
    else
        INDEL_BAN_VAL=$INDEL_BAN_ARG # for later use
        if [[ $INDEL_BAN_ARG =~ ban ]] ; then
            INDEL_BAN_ARG="[DI]"
        elif [[ $INDEL_BAN_ARG -lt 10 ]] ; then
            INDEL_BAN_ARG=$(echo "[1-$INDEL_BAN_ARG][DI]")
        elif [[ $INDEL_BAN_ARG -lt 20 ]] ; then
            INDEL_BAN_ARG=$(echo "[^1-9][1-9][DI]|^[1-9][DI]|[1][0-$((INDEL_BAN_ARG-10))][DI]")
        elif [[ $INDEL_BAN_ARG -lt 30 ]] ; then
            INDEL_BAN_ARG=$(echo "[^1-9][1-9][DI]|^[1-9][DI]|[1-2][0-$((INDEL_BAN_ARG-20))][DI]")
        else
            INDEL_BAN_ARG=$(echo "[^1-9][1-9][DI]|^[1-9][DI]|[1-2][0-9][DI]|30[DI]")
        fi
    fi
fi
if [[ $BAM_LIST == " " ]] ; then
     echo "No BAM file(s) provided!"
     exit_err
else
    for BAM in $BAM_LIST
    do
        if ! [[ $BAM =~ ^.+\.bam$ || $BAM =~ ^.+\.sam$ ]] ; then # add verification by samtools
            echo $BAM "is not a valid input file name!"
            exit_err
        fi
    done
    samtools view -@ $((THR_ARG-1)) $(echo $BAM_LIST | cut -d" " -f1) | head -100 | gawk '{ if (and($2, 0x1)) { exit 0 } else { exit 1 } }'
    if [ $? -eq 0 ]; then
        PAIRED_ARG=1 # paired-end
    elif [ $? -eq 1 ]; then
        PAIRED_ARG=0 # single-end
    fi
fi

if [ $? -eq 0 ]; then
    PAIRED_ARG=1 # paired-end
elif [ $? -eq 1 ]; then
    PAIRED_ARG=0 # single-end
fi

if [[ $MODE_ARG =~ ^.*slamdunk.*$ && $PAIRED_ARG == 1 ]]; then
    echo "Currently \"slamdunk\" mode requires single-end input and provided BAM file(s) contain(s) paired-end sequences. Rebound will attempt to convert to single-end retaining original mapping positions. This should work, but is rather experimental."
fi

if ! [ -f "$REF_ARG" ] ; then
    echo "Fatal error: $REF_ARG not found."
    exit_err
elif ! [[ $REF_ARG =~ ^.+\.fasta$ || $REF_ARG =~ ^.+\.fa$ ]] ; then
    echo "$REF_ARG does not look like a valid reference file."
    exit_err
else
    echo "Genome reference: $REF_ARG"
fi

if ! [ -f "$GTF_ARG" ] ; then
    echo "Fatal error: $GTF_ARG not found."
    exit_err
elif ! [[ $GTF_ARG =~ ^.+\.gtf$ || $REF_ARG =~ ^.+\.gff$ ]] ; then
    echo "$GTF_ARG does not look like a valid annotation file."
    exit_err
else
    echo "Annotations: $GTF_ARG"
fi

if [[ $NODUPS_ARG == 1 ]]; then
    echo "Duplicate removal [ -d ] is ON"
fi

if [[ $STRANDEDNESS_ARG =~ reverse|forward|unstranded ]] ; then # change if statemet to case statement
    echo "Strandedness is set to $STRANDEDNESS_ARG."
    if [[ $STRANDEDNESS_ARG =~ forward ]] ; then
        STRANDEDNESS_ARG="yes"
    elif [[ $STRANDEDNESS_ARG =~ reverse ]] ; then
        STRANDEDNESS_ARG="reverse"
    elif [[ $STRANDEDNESS_ARG =~ unstranded ]] ; then
        STRANDEDNESS_ARG="no"
    fi
else
    echo "Fatal error: Strandedness must be one of \"forward\", \"reverse\" or \"unstranded\"."
    exit_err
fi

if ! [[ $MODE_ARG =~ ^.*htseq_only.*$ ]] ; then
    if ! [[ -z $QUALITY_ARG ]] ; then
        MAPQ_FILTER=$(echo $QUALITY_ARG | cut -d"," -f1)
        NM_FILTER=$(echo $QUALITY_ARG | cut -d"," -f2)
        MISMATCH_FILTER=$(echo $QUALITY_ARG | cut -d"," -f3)
        MISMATCH_QUALITY=$(echo $QUALITY_ARG | cut -d"," -f4)
        MISMATCH_QUALITY_TC=$(echo $QUALITY_ARG | cut -d"," -f5)
    elif [[ $MAPPER_ARG =~ ^.*star.*$ ]] ; then # defaults star
        MAPQ_FILTER=254
        NM_FILTER=20
        MISMATCH_FILTER="0.95"
        MISMATCH_QUALITY=27
        MISMATCH_QUALITY_TC=27
    elif [[ $MAPPER_ARG =~ ^.*bbmap.*$ ]] ; then # defaults bbmap
        MAPQ_FILTER=20
        NM_FILTER=100000 # off
        MISMATCH_FILTER="0.95"
        MISMATCH_QUALITY=27
        MISMATCH_QUALITY_TC=27
    fi
    printf " ------------------------------------\n Read filter settings\n Min mapping quality (MAPQ): ${BOLD}${MAPQ_FILTER}${NORM}\n Max edit distance (NM): ${BOLD}${NM_FILTER}${NORM}\n Min fraction of matching bases (XF): ${BOLD}${MISMATCH_FILTER}${NORM}\n Min mismatch call quality (QV): ${BOLD}${MISMATCH_QUALITY}${NORM}\n Min mismatch call quality TC (QVTC): ${BOLD}${MISMATCH_QUALITY_TC}${NORM}\n"
    # echo $INDEL_BAN_VAL
    if ! [[ -z $INDEL_BAN_ARG ]] ; then
        if [[ $INDEL_BAN_VAL =~ disabled ]] ; then
            printf " ----------------------------------\n"
        elif [[ $INDEL_BAN_VAL =~ ban ]] ; then
            printf " Ignoring reads with indels.\n ------------------------------------\n"
        else
            printf " Ignoring reads with indels ${BOLD}<= ${INDEL_BAN_VAL}${NORM} nucleotides.\n ------------------------------------\n"
        fi
    else
        printf " ----------------------------------\n"
    fi

    XI_VALUE=$(echo ${MISMATCH_FILTER}*100|bc)
    REP1="Min mapping quality (MAPQ):\t${MAPQ_FILTER}"
    REP2="Max edit distance (NM):\t${NM_FILTER}"
    REP3="Min fraction of matching bases (XF):\t${XI_VALUE}%%"
    REP4="Min mismatch call quality (QV):\t${MISMATCH_QUALITY}"
    REP5="Min mismatch call quality TC (QVTC):\t${MISMATCH_QUALITY_TC}"
    REP6="Indel ban threshold:\t${INDEL_BAN_VAL}"

fi

if ! [[ $THR_ARG =~ ^[0-9]+$ ]] ; then
    echo "Threads number must be an integer!"
    exit_err
fi

echo "CPU threads: $THR_ARG"
echo "RAM per thread: $MEM_ARG"

##### main routine #####

START_TIME=$SECONDS

if [[ $MODE_ARG == "htseq_only" ]]; then # test mode
    for BAM in $BAM_LIST
    do
        BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' )
        BAM_FILE=$( echo $BAM | xargs basename)
            echo "Running only HTSseq. Irrelevant parameters will be ignored."
            htseq-count -o ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.sam -r pos -a 0 -f bam -s $([[ $STRANDEDNESS_ARG == "forward" ]] && echo "yes" || echo "reverse") --nonunique $([[ $UNIQ_ARG -eq 1 ]] && echo "none" || echo "all") -i gene_id $BAM $GTF_ARG > ${OUTDIR_ARG}/${BAM_NAME}_counts.txt
            samtools view -H $BAM > ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam
            cat ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.sam >> ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam
            samtools view -@ $((THR_ARG-1)) -hb ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam > ${OUTDIR_ARG}/${BAM_NAME}_hts_only.bam && rm ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam
            echo "Finished!"
    done
    echo "All completed. Exiting."
    exit 0
fi

for BAM in $BAM_LIST
do
    BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' )
    BAM_FILE=$( echo $BAM | xargs basename )
#####################################
    printf " -----------------------------------\n Processing%s $BAM_FILE\n"

    if [[ $PAIRED_ARG == 1 ]]; then
        echo " Data is paired-end"
        CLIP_ARG=1
    else
        echo " Data is single-end"
        CLIP_ARG=0
    fi

    printf " -----------------------------------\n"
#####################################
    if [[ $MD_ARG -eq 0 ]] ; then
        samtools view -@$((THR_ARG-1)) $BAM | head -1 | gawk '{if ($0 ~ /^.*MD:Z:.*$/){exit 0}}'
        if ! [ $? -eq 0 ]; then
            echo "MD tags not found and will be added with \"samtools calmd\"."
            MD_ARG=1
        fi
    else # force MD
        MD_ARG=1
    fi
#####################################
    if [[ $MODE_ARG == "normal" ]] ; then
        echo "Annotating and counting reads."
        htseq-count -o ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.sam -r pos -a 0 -f bam -s $STRANDEDNESS_ARG --nonunique $([[ $UNIQ_ARG -eq 1 ]] && echo "none" || echo "all") -i gene_id $BAM $GTF_ARG > $OUTDIR_ARG/${BAM_NAME}_counts.txt
        samtools view -H $BAM > ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam
        if [[ $MAPPER_ARG == "bbmap" ]] ; then
            cat ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.sam | gawk 'BEGIN{ FS="\t"; OFS="\t"} ; { if ($0 ~ /^.*YI:j:.*$/ && $6 ~ /M/){
            gsub("YI:j:", "YI:f:", $0); print $0 }}' >> ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam # remove or modify later, option specific for bbmap; this also checks if CIGARs exist
        elif [[ $MAPPER_ARG == "star" ]] ; then
            cat ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.sam | gawk 'BEGIN{ FS="\t"; OFS="\t"} ; { if ($6 ~ /M/) print $0 }' >> ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam # remove or modify later option specific for bbmap; this also checks CIGARs
        fi
        BAM="${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam"
    elif [[ $MODE_ARG == "slamdunk" ]] ; then
        if ! [[ -z $STRANDEDNESS_ARG ]] ; then
            echo "[ -s ] flag will be ignored in mode \"slamdunk\"."
        fi
        if ! [[ -z $UNIQ_ARG ]] ; then
            echo "[ -u ] flag will be ignored in mode \"slamdunk\"."
        fi
    elif [[ $MODE_ARG == "nocount" ]] ; then # test mode
        echo "Skipping counting with HTseq."
    fi
#####################################
    rm -f $OUTDIR_ARG/${BAM_NAME}_hts_strand.sam
#####################################
    samtools view -@$((THR_ARG-1)) $BAM | head -10000 | gawk '{if ($6 ~ /N/){exit 1}}'

    if [[ $? -eq 1 ]]; then
        echo "i(N)trons in CIGARs will be re-labeled as (D)eletions"
        NINCIG=1
    fi
    echo "Preparing SLAM-seq tags"
        # only including mapped reads (0x4 flag inverse)
        # echo $INDEL_BAN_ARG
        # echo $NODUPS_ARG
    if [[ -z $INDEL_BAN_ARG || $INDEL_BAN_VAL =~ disabled ]] && [[ $NODUPS_ARG == 1 ]] ; then
        samtools view -@$((THR_ARG-1)) -h -F 0x4 -F 0x400 $BAM 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log
    elif [[ -z $INDEL_BAN_ARG || $INDEL_BAN_VAL =~ disabled ]] && [[ $NODUPS_ARG == 0 ]] ; then
        samtools view -@$((THR_ARG-1)) -h -F 0x4 $BAM 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log
    elif [[ ! -z $INDEL_BAN_ARG && $NODUPS_ARG == 1 ]] ; then
        samtools view -@$((THR_ARG-1)) -h -F 0x4 -F 0x400 $BAM 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log | gawk -v INDEL_BAN_ARG=$INDEL_BAN_ARG -v BAM_NAME=$BAM_NAME -v OUTDIR_ARG=$OUTDIR_ARG 'BEGIN { FS="\t"; OFS="\t" } ; {if ($1 ~ /^@/) { print $0 } else { stats["total"]++ ; if ($6 !~ INDEL_BAN_ARG) { stats["indel_filter_pass"]++ ; print $0 } }} ; END { print "Total reads (excl. dups and unmapped):", stats["total"] > OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log" ; print "Reads passing indel filter:", stats["indel_filter_pass"] > OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log" }'
    elif [[ ! -z $INDEL_BAN_ARG && $NODUPS_ARG == 0 ]] ; then
        samtools view -@$((THR_ARG-1)) -h -F 0x4 $BAM 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log | gawk -v INDEL_BAN_ARG=$INDEL_BAN_ARG -v BAM_NAME=$BAM_NAME -v OUTDIR_ARG=$OUTDIR_ARG 'BEGIN { FS="\t"; OFS="\t" } ; {if ($1 ~ /^@/) { print $0 } else { stats["total"]++ ; if ($6 !~ INDEL_BAN_ARG) { stats["indel_filter_pass"]++ ; print $0 } }} ; END { print "Total reads (excl. unmapped):", stats["total"] > OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log" ; print "Reads passing indel filter:", stats["indel_filter_pass"] > OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log" }'
    fi |\

    if [[ $MD_ARG -eq 1 ]]; then
        # if [[ $(samtools view -H $BAM | head -1) =~ ^.*oordinate.*$ ]]; then
        #     cat
        # else
        #     samtools sort -m $MEM_ARG -@$((THR_ARG-1)) -T $OUTDIR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log
        # fi |\
        samtools sort -m $MEM_ARG -@$((THR_ARG-1)) -T $OUTDIR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log |\
        samtools calmd -@ $((THR_ARG-1)) - $REF_ARG 2> /dev/null
        # 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr_tmp.log # overwrite to save space
    else
        cat
    fi |\

    if [[ $CLIP_ARG == 1 ]]; then
        samtools sort -n -m $MEM_ARG -@$((THR_ARG-1)) -T $OUTDIR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log |\
        bam clipOverlap --in - --out - --readName --storeOrig OC --stats 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_clipOverlap_stderr.log |\
        samtools sort -m $MEM_ARG -@$((THR_ARG-1)) -T $OUTDIR_ARG -O sam - 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log
    else
        cat
    fi |\

    if [[ $NINCIG == 1 ]]; then
        gawk 'BEGIN{FS="\t"; OFS="\t"; print "Replacing i(N)trons by (D)eletions in CIGARs" > "/dev/stderr"}; function tag_finder(tag, i){for (i=1; i<=NF+1; i++){if (i==NF+1) {print "Fatal error: unable to find " tag " at record " $0; exit 1 } else if (substr($i,1,5)!=tag) {continue} else { return $i }}}; {if ($1 ~ /^@/) { print $0 } else { gsub($6, gensub("N","D","g",$6), $0); if ($0 ~ /^.*NM:i:.*$/) { print $0, "ED:i:" gensub("NM:i:","","g",tag_finder("NM:i:")) } else { print $0 }}}' |\
        samtools calmd -@ $((THR_ARG-1)) - $REF_ARG 2> /dev/null
        # 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr_tmp.log # overwrite to save space
    else
        if [[ $MODE_ARG == "slamdunk" && $MD_ARG -eq  0 ]] ; then
            cat # the only rare case MD-tag addition is skipped
        else
            samtools calmd -@ $((THR_ARG-1)) - $REF_ARG 2> /dev/null
            # 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr_tmp.log # overwrite to save space
        fi
    fi |\

    if [[ $MODE_ARG == "slamdunk" && $PAIRED_ARG == 1 ]]; then
        gawk 'BEGIN{FS="\t"; OFS="\t"; print "Converting to single-end BAM" > "/dev/stderr" }; {if ($1 ~ /^@/) { print $0 } else {if ($7 != "*") {$7 = "*"; $8 = 0; $9 = 0; if ((and($2, 0x20) && and($2, 0x80)) || (and($2, 0x10) && and($2, 0x40))) {$2 = 0} else if ((and($2, 0x10) && and($2, 0x80)) || (and($2, 0x20) && and($2, 0x40))) {$2 = 16} print $0 }}}'
    else
        cat
    fi |\
    # tee $OUTDIR_ARG/${BAM_NAME}_prior_to_rebound.sam |\
    gawk -v NODUPS_ARG=$NODUPS_ARG -v MAPPER_ARG=$MAPPER_ARG -v MODE_ARG=$MODE_ARG -v MAPQ_FILTER=$MAPQ_FILTER -v NM_FILTER=$NM_FILTER -v MISMATCH_FILTER=$MISMATCH_FILTER -v MISMATCH_QUALITY_TC=$MISMATCH_QUALITY_TC -v MISMATCH_QUALITY=$MISMATCH_QUALITY -v REBOUND_COMMAND="$REBOUND_COMMAND" -v BAM_NAME="$BAM_NAME" -v OUTDIR_ARG="$OUTDIR_ARG" -v STRANDEDNESS_ARG="$STRANDEDNESS_ARG" -f $AWK_FILE 2> $OUTDIR_ARG/stderr/${BAM_NAME}_rebound_stderr.log 1> $OUTDIR_ARG/${BAM_NAME}_rebound.sam
    REBOUND_EXIT=$?

    # head -1000 $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr_tmp.log > $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr.log && rm -f $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr_tmp.log

    # if [[ $(cat $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr.log | wc -l) -eq 1000 ]]; then
    #    echo "File truncated here due to size !!! ######### !!! ######### !!!" >>  $OUTDIR_ARG/stderr/${BAM_NAME}_calmd_stderr.log
    #    echo "${BAM_NAME}_calmd_stderr.log file was truncated due to size"
    # fi

    # rm -f $OUTDIR_ARG/stderr/calmd_stderr_tmp.log
    rm -f ${OUTDIR_ARG}/${BAM_NAME}_tmp_htseq.sam
    # samtools view -@ $((THR_ARG-1)) -hb $OUTDIR_ARG/${BAM_NAME}_prior_to_rebound.sam > $OUTDIR_ARG/${BAM_NAME}_prior_to_rebound.bam && rm $OUTDIR_ARG/${BAM_NAME}_prior_to_rebound.sam

    samtools view -@ $((THR_ARG-1)) -hb ${OUTDIR_ARG}/${BAM_NAME}_rebound.sam > ${OUTDIR_ARG}/${BAM_NAME}_rebound.bam && rm ${OUTDIR_ARG}/${BAM_NAME}_rebound.sam
    samtools index ${OUTDIR_ARG}/${BAM_NAME}_rebound.bam
    if [[ $REBOUND_EXIT -eq 0 ]]; then
        echo "Processing $BAM_FILE successful"
        echo "Output BAM saved to $OUTDIR_ARG/${BAM_NAME}_rebound.bam"
        if ! [[ $MODE_ARG =~ ^.*htseq_only.*$ ]] ; then
            printf "${REP1}\n" >> ${OUTDIR_ARG}/${BAM_NAME}_rebound_stats.log
            printf "${REP2}\n" >> ${OUTDIR_ARG}/${BAM_NAME}_rebound_stats.log
            printf "${REP3}\n" >> ${OUTDIR_ARG}/${BAM_NAME}_rebound_stats.log
            printf "${REP4}\n" >> ${OUTDIR_ARG}/${BAM_NAME}_rebound_stats.log
            printf "${REP5}\n" >> ${OUTDIR_ARG}/${BAM_NAME}_rebound_stats.log
            printf "${REP6}\n" >> ${OUTDIR_ARG}/${BAM_NAME}_rebound_stats.log
        fi
    else
        echo "${BOLD}Processing $BAM_FILE FAILED${NORM}"
    fi
done

TOTAL_TIME=$(($SECONDS - $START_TIME))
echo "Finished processing" $(echo "$BAM_LIST" | wc -w | awk '{print $1}') "file(s) in" $(echo "scale=2;$TOTAL_TIME/60" | bc) "minute(s)."
