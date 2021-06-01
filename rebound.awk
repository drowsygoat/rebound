BEGIN{
start_time = systime()
OFS="\t"
FS="\t"
for (n=0;n<256;n++){
    phred_conv[sprintf("%c",n+33)]=n
}
print "S(L)AM converter started @ " strftime() > "/dev/stderr"
if (MODE_ARG == "normal"){
    if (NODUPS_ARG == 1){
        print "gene_id", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "DUP" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
    }else{
        print "gene_id", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
    }
}

# --- FUNCTIONS --- #

function tag_finder(tag,    i){
    if (substr($0,1,1)!="@"){
        for (i=1; i<=NF+1; i++){
            if (i==NF+1) {
                print "Fatal error: unable to find " tag " at record " $0 > "/dev/stderr"
                exit 1
            }else if (substr($i,1,5)!=tag){
                continue
            }else{
                return $i
            }
        }
    }
}

function rmcol(col,     i){
    for (i=col; i<NF; i++){
        $i = $(i+1)
    }
    NF--
}

function MD_extender(md,    out,a,b) {
    out=""
    match(md,/MD:Z:/)
    a=substr(md,RSTART+RLENGTH)
    while(match(a,/\^[A-Z]+/)) {
        a=substr(a,1,RSTART-1) "h" substr(a,RSTART+RLENGTH)
    }
    n = patsplit(a, b, /[0-9]+/, seps)
    for (i=1;i<=n;i++) {
        for(j=1;j<=b[i];j++) {
            out=out "."
        }
        if (seps[i] ~ /[ACTGN]/) {
            out=out seps[i]
        }else if (seps[i] ~ /h/) {
                out=out
        }
    }
    return out
}
# returns matched + mismatched (M) portion of the reads and corresponding quality string
function MM_extractor(read, quali, ec,     a,b,j,i,k) {
    mm_read=""; mm_quali=""
    split(quali, c, "")
    split(read, b, "")
    split(ec, a, "")
    j=1
    for (i in a) {
        if (a[i] ~ /[I]/) {
            j++
            continue
        }else if (a[i] ~ /[S]/) {
            j++
            continue
        }else if (a[i] ~ /[DN]/) {
            continue
        }else if (a[i] ~ /[M]/) {
            mm_read = mm_read b[j]
            mm_quali = mm_quali c[j]
            j++
            continue
        }else{
            print "Fatal error: Unexpected character(s) in CIGAR string of read" $0 > "/dev/stderr"
            exit 1
        }
    }
}

function cigar_extender(cygaro,   out, b,i,j){
    # expand CIGAR string into a longer string (out) with a single letter per alignment position; works for CIGARs with [SIDNM] operations
    out = ""
    n = patsplit(cygaro, b, /[^[:digit:]]/, seps)
    # loop through CIGAR "segments" to get their lengths (--> seps)
    for (i=0; i<=n-1; i++) {
        # iteratively add letter codes correspondng to CIGAR operations of CIGAR string being expanded
       for (j=1; j<=seps[i]; j++) {
           out=out b[i+1]
       }
   }
   return out
}

function MP_producer(MD_input, ec, read,      a, b){
    # outputs array (read_mis) with mismatched bases in the read, and their respective indices; if there are not gaps in the alignment, and if nothing was clipped from the read, these indices are identical to the indices of the mismatched reference bases; indexing starts from the start of the alignment;
    # arrays ref_mis and read_mis contain mismatched bases and their positions (with respect to start of the alignment); an offset array is populated that contains offset values for consecutive pairs; in case of reads that have no I, D, N, or S in their CIGARs the values in read_mis and ref_mis are identical
    sum_md=0; delete(ref_mis); delete(read_mis)
    match(MD_input,/MD:Z:/)
    a=substr(MD_input,RSTART+RLENGTH)
    while(match(a,/\^[A-Z]+/)) {
        a=substr(a,1,RSTART-1) "h" RLENGTH-1 "^" substr(a,RSTART+RLENGTH)
    }
    n=patsplit(a, b, /[[:alpha:]^]/, seps)
    for(i=1; i<=n; i++){
        # print "sepsi", seps[i], i
        if (seps[i] == ""){
            seps[i] = 0
        }
        if (b[i]=="^"){
            sum_md+=seps[i-1]
        }else if (b[i]=="h"){
            sum_md+=seps[i-1]
        }else{
            sum_md+=seps[i-1]+1
        }
        if (b[i]!="^" && b[i]!="h"){
            ref_mis[sum_md]=b[i]
        }
    }
    for (i in ref_mis){
        j=0
        k=0
        for (m=1; m<=length(ec); m++){
            if(substr(ec,m,1) == "M"){
                k++
                if (k==i*1){
                    break
                }
            }else if (substr(ec,m,1) == "S" || substr(ec,m,1) == "I"){
                j++
            }else if (substr(ec,m,1) == "D" || substr(ec,m,1) == "N"){
                j--
                k++ # only deletions/introns need to be acounted for to match the positions in CIGARs
            }
        }

        offset[i] = j
        read_mis[i+j] = substr(read,i+j,1)

        mp = "MP:Z:"
        sep = ","
        cnt=1 # counter
        for (i in ref_mis) {
            j = offset[i]
            # if (ref_mis[i] == "A" && read_mis[i+j] == "A"){
            #     mp = mp "0:" i+j ":" i
            # }
            if (ref_mis[i] == "A" && read_mis[i+j] == "C"){
                mp = mp "1:" i+j ":" i
            }
            else if (ref_mis[i] == "A" && read_mis[i+j] == "G"){
                mp = mp "2:" i+j ":" i
            }
            else if (ref_mis[i] == "A" && read_mis[i+j] == "T"){
                mp = mp "3:" i+j ":" i
            }
            else if (ref_mis[i] == "A" && read_mis[i+j] == "N"){
                mp = mp "4:" i+j ":" i
            }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "A"){
                mp = mp "5:" i+j ":" i
            }
            # else if (ref_mis[i] == "C" && read_mis[i+j] == "C"){
            #     mp = mp "6:" i+j ":" i
            # }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "G"){
                mp = mp "7:" i+j ":" i
            }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "T"){
                mp = mp "8:" i+j ":" i
            }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "N"){
                mp = mp "9:" i+j ":" i
            }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "A"){
                mp = mp "10:" i+j ":" i
            }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "C"){
                mp = mp "11:" i+j ":" i
            }
            # else if (ref_mis[i] == "G" && read_mis[i+j] == "G"){
            #     mp = mp "12:" i+j ":" i
            # }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "T"){
                mp = mp "13:" i+j ":" i
            }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "N"){
                mp = mp "14:" i+j ":" i
            }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "A"){
                mp = mp "15:" i+j ":" i
            }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "C"){
                mp = mp "16:" i+j ":" i
            }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "G"){
                mp = mp "17:" i+j ":" i
            }
            # else if (ref_mis[i] == "T" && read_mis[i+j] == "T"){
            #     mp = mp "18:" i+j ":" i
            # }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "N"){
                mp = mp "19:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "A"){
                mp = mp "20:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "C"){
                mp = mp "21:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "G"){
                mp = mp "22:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "T"){
                mp = mp "23:" i+j ":" i
            }
            # else if (ref_mis[i] == "N" && read_mis[i+j] == "N"){
            #     mp = mp "24:" i+j ":" i
            # }
            if (cnt < length(ref_mis)){
                mp = mp ","
                cnt++
            }
        }
        if (length(ref_mis) == 0){
            print "This should not have happened" > "/dev/stderr"
            exit 1
        }
    }
}

function TCRA_producer(mm_read_, mm_quali_, md_extended,    i){
    delete count; delete mmms; delete quali
    split(mm_read_, mm_read_split, "")
    split(mm_quali_, mm_quali_split, "")
    split(md_extended, md_split, "")
    if (length(mm_read_split) != length(mm_quali_split)){
        print "Strings of different lengths. This should not have happened. Record:" > "/dev/stderr"
        print $0 > "/dev/stderr"
        exit 1
    }
    for (i=1; i<=length(mm_read_split); i++) {
        if (md_split[i] == "."){
            mmms[i] = mm_read_split[i] mm_read_split[i]
        }else{
            mmms[i] = md_split[i] mm_read_split[i]
        }
    }
    if (length(mmms) != length(mm_quali_split)){
        print "Strings of different lengths. This should not have happened. Record:" > "/dev/stderr"
        print $0 > "/dev/stderr"
        exit 1
    }
    for (i in mmms) {
        if (phred_conv[mm_quali_split[i]] >= MISMATCH_QUALITY){
            ###### A -> X
            }if (mmms[i] == "AG"){
                if ((and($2, 0x10) && and($2, 0x80)) || (and($2, 0x20) && and($2, 0x40)) || $2 == 16){
                    if (phred_conv[mm_quali_split[i]] >= MISMATCH_QUALITY_TC){
                        count["AG"]++
                    }
                }else{
                    count["AG"]++
                }
            }else if (mmms[i] == "AA"){
                count["AA"]++
            }else if (mmms[i] == "AC"){
                count["AC"]++
            }else if (mmms[i] == "AT"){
                count["AT"]++
            }else if (mmms[i] == "AN"){
                count["AN"]++
                # quali["AN"]=quali["AN"]+(1-10^(-(phred_conv[mm_quali_split[i]]/10)))

            ###### C -> X
            }else if (mmms[i] == "CA"){
                count["CA"]++
            }else if (mmms[i] == "CC"){
                count["CC"]++
            }else if (mmms[i] == "CG"){
                count["CG"]++
            }else if (mmms[i] == "CT"){
                count["CT"]++
            }else if (mmms[i] == "CN"){
                count["CN"]++
                # quali["CN"]=quali["CN"]+(1-10^(-(phred_conv[mm_quali_split[i]]/10)))

            ###### G -> X
            }else if (mmms[i] == "GA"){
                count["GA"]++
            }else if (mmms[i] == "GC"){
                count["GC"]++
            }else if (mmms[i] == "GG"){
                count["GG"]++
            }else if (mmms[i] == "GT"){
                count["GT"]++
            }else if (mmms[i] == "GN"){
                count["GN"]++
                # quali["GN"]=quali["GN"]+(1-10^(-(phred_conv[mm_quali_split[i]]/10)))

            ###### T -> X
            }else if (mmms[i] == "TA"){
                count["TA"]++
            }else if (mmms[i] == "TC"){
                if ((and($2, 0x20) && and($2, 0x80)) || (and($2, 0x10) && and($2, 0x40)) || $2 == 0){
                    if (phred_conv[mm_quali_split[i]] >= MISMATCH_QUALITY_TC){
                        count["TC"]++
                    }
                }else{
                    count["TC"]++
                }
            }else if (mmms[i] == "TG"){
                count["TG"]++
            }else if (mmms[i] == "TT"){
                count["TT"]++
            }else if (mmms[i] == "TN"){
                count["TN"]++
                # quali["TT"]=quali["TT"]+(1-10^(-(phred_conv[mm_quali_split[i]]/10)))

            ###### N -> X
            }else if (mmms[i] == "NA"){
                count["NA"]++
            }else if (mmms[i] == "NC"){
                count["NC"]++
            }else if (mmms[i] == "NG"){
                count["NG"]++
            }else if (mmms[i] == "NT"){
                count["NT"]++
            }else if (mmms[i] == "NN"){
                count["NN"]++
            # print "true"
        # }else{
            # print "false"
        }
    }
    # write TC tags, only needed for slamdunk mode
    # reverse strand + first in pair or forward strand + second in pair or forward unpaired
    if ((and($2, 0x20) && and($2, 0x80)) || (and($2, 0x10) && and($2, 0x40)) || $2 == 0){
        tc = count["TC"]+0
    #reverse strand + second in pair or forward strand + first in pair or reverse unpaired
    }else if ((and($2, 0x10) && and($2, 0x80)) || (and($2, 0x20) && and($2, 0x40)) || $2 == 16){
        tc = count["AG"]+0
    }
    sp = ","
    ra =      count["AA"]+0 sp count["AC"]+0 sp count["AG"]+0 sp count["AT"]+0 sp count["AN"]+0 sp count["CA"]+0 sp count["CC"]+0 sp count["CG"]+0 sp count["CT"]+0 sp count["CN"]+0 sp count["GA"]+0 sp count["GC"]+0 sp count["GG"]+0 sp count["GT"]+0 sp count["GN"]+0 sp count["TA"]+0 sp count["TC"]+0 sp count["TG"]+0 sp count["TT"]+0 sp count["TN"]+0 sp count["NA"]+0 sp count["NC"]+0 sp count["NG"]+0 sp count["NT"]+0 sp count["NN"]+0 sp

    if (count["AC"] + count["AG"] + count["AT"] + count["AN"] + count["CA"] + count["CG"] + count["CT"] + count["CN"] + count["GA"] + count["GC"] + count["GT"] + count["GN"] + count["TA"] + count["TC"] + count["TG"] + count["TN"] + count["NA"] + count["NC"] + count["NG"] + count["NT"] + count["NN"] > 0){

        xi = (count["AA"] + count["CC"] + count["GG"] + count["TT"] + 0) / (count["AA"] + count["AC"] + count["AG"] + count["AT"] + count["AN"] + count["CA"] + count["CC"] + count["CG"] + count["CT"] + count["CN"] + count["GA"] + count["GC"] + count["GG"] + count["GT"] + count["GN"] + count["TA"] + count["TC"] + count["TG"] + count["TT"] + count["TN"] + count["NA"] + count["NC"] + count["NG"] + count["NT"] + count["NN"])
        # keeping Ns, as they may signify a bad read (?)
    }else{ # there reads with zero counts will not be included anyways
        xi = 0
    }
}

function TCRA_producer_simple(read, ec,     a, count) {
    # use this for reads with no mismatches and also if reads were entirely clipped (cases of 100% overlap) == only S in CIGARS; this reads may actually be tossed completely, but they can be used for DE analysis when one wants to use exactly the same pool of reads for both, DE and SLAM-seq; in such sase, the reads would need to get the identi score of 1, otherwise will be filtered out; currenyl these reads are included in the resulting BAM file, but not in the counts, and they don't contribute to T coverage
    if ($6 ~ /^[^M]+$/){ # no matches and mismatches == full clip
        print "Warning: first condition in TC simple function is met at record: " $0 > "/dev/stderr"
        tc_simple = 0
        ra_simple = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
        xi_simple = 0 # may be set to 1 if filtering of reads with no Ms in CIGARs is not desired
    }else{
        # CIGAR case with M only # switch if statement order here
        if ($6 ~ /^[^SDNI]+$/) {
            a = $10
        # CIGAR cases with I,N,D,S
        }else{
            MM_extractor($10, $11, ec)
            a = mm_read
        }
        delete(count)
        split(a, a_split, "")
        for (i in a_split){
            if       (a_split[i] == "A"){
                count["A"]++
            }else if (a_split[i] == "C"){
                count["C"]++
            }else if (a_split[i] == "G"){
                count["G"]++
            }else if (a_split[i] == "T"){
                count["T"]++
            }else if (a_split[i] == "N"){
                count["N"]++
            }
        }
        tc_simple = 0
        sp = ","
        ra_simple = count["A"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["C"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["G"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["T"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["N"]+0 sp
        ra_rc_simple = count["T"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["G"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["C"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["A"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["N"]+0 sp
        xi_simple = 1
    }
}

# ---- main routine ---- #

{
    if ($1 ~ /^@/){
        header=1
        if ($1 ~ /^@RG/){
            suf="DS:{'sequenced':0,'mapped':0,'filtered':0,'mqfiltered':0,'idfiltered':0,'nmfiltered':0,'multimapper':0,'dedup':0,'snps':0,'annotation':'','annotationmd5':''}" # irrelevant, update to include relevant data from sample table
            # rg_line=$1 ", ID:" BAM_NAME ", SM:" BAM_NAME ":NA:NA," suf
            print $1,"ID:"BAM_NAME,"SM:"BAM_NAME":NA:NA", suf
        }else{
            print $0
        }
    }
    if ($1 !~ /^@/){
        ############################## Determine strand ###################################
        if ((and($2, 0x20) && and($2, 0x80)) || (and($2, 0x10) && and($2, 0x40)) || $2 == 0){
            strand2 = "forward"
        }else if ((and($2, 0x10) && and($2, 0x80)) || (and($2, 0x20) && and($2, 0x40)) || $2 == 16){
            strand2 = "reverse"
        }else{
            strand2 = "undetermined"
        }
        ############################## Get tags ###################################
        if ($0 ~ /^.*MD:Z:.*$/){ # get MD tag field if present
            md = tag_finder("MD:Z:")
            md_value = gensub(/MD:Z:/, "", "g", md)
        }
        if ($0 ~ /^.*XF:Z:.*$/){ # get XF tag field if present
            xf = tag_finder("XF:Z:")
            xf_value = gensub(/XF:Z:/, "", "g", xf)
        }
        if (MAPPER_ARG == "star"){
            if ($0 ~ /^.*ED:i:.*$/){ # get ED tag field if present
                ed = tag_finder("ED:i:")
                ed_value = gensub(/ED:i:/, "", "g", ed)
            }
            ############################## Convert tags ###################################
            for (i=1; i<=NF+1; i++){ # converting ED tag to recover NM tag (original NM tags Ns were replaced with Ds, and recalculated NM values were inflated)
                if (substr($i,1,5)=="NM:i:"){
                    $i=ed
                    break
                }else if (i==NF+1) {
                    print "Error: unable to find NM tag at record " $0 > "/dev/stderr"
                    exit 1
                }
            }
            for (i=1; i<=NF+1; i++){ # removing ED tag after converting it to NM tag
                if (substr($i,1,5)=="ED:i:"){
                    rmcol(i)
                    break
                }else if (i==NF+1) {
                    print "Error: unable to find ED tag at record " $0 > "/dev/stderr"
                    exit 1
                }
            }
        }else if (MAPPER_ARG == "bbmap"){
            if ($0 ~ /^.*YI:f:.*$/){ # get YI tag field if present
                ed = tag_finder("YI:f:")
                ed_value = gensub(/YI:f:/, "", "g", ed)
                ed_value = ed_value/100
            }
        }
    }
    # skip header and modify @RG by adding DS. SM can be supplied as argument.
    if ($1 !~ /^@/ && $6 !~ /^[^M]+$/){ # only include reads with  matches to the template (many reads like this are a result of overlap clipping)
        stats["reads_with_M_in_cig"]++
        if (header==1){ # first non-header record
            header=0
            print "@PG","ID:rebound.sh","VN:1.0","CL:"REBOUND_COMMAND
        }
        # if (and($2, 0x4)){ # commented to exclude unmapped
        #     print $0,"XI:f:0"
        if (md_value ~ /^[0-9]*(\^[ACTGN]+[0-9]+)*$/){
            stats["reads_without_mismatches"]++
            ec=cigar_extender($6)
            TCRA_producer_simple($10, ec)
            # if (MAPPER_ARG == "bbmap"){
            #     ed_value=xi_simple
            # }
            if ($5 >= MAPQ_FILTER && ed_value <= NM_FILTER && xi_simple >= MISMATCH_FILTER){
                stats["reads_without_mismatches_filtered"]++
                if (MODE_ARG == "slamdunk"){
                    print $0,"XI:f:"xi_simple,"TC:i:"tc_simple,"RA:Z:"ra_simple
                }else if (MODE_ARG == "normal"){
                    print $0,"XI:f:"xi_simple,"TC:i:"tc_simple,"RA:Z:"ra_simple
                    if (xf !~ /_ambiguous|_no_feature|_too_low_aQual|_not_aligned|_alignment_not_unique/){  # move this filter before tha main loop to avoid processing of reads that will be excluded anyways by HTseq
                        split(xf_value,a,"__") # only if annotation is appended with strandedness
                        gene_id = a[1]; strand = a[2]
                        split(ra_simple,b,",")
                        if ((strand == "+" && strand2 == "reverse") || (strand == "-" && strand2 == "forward")){ # make a conditional to print only first record like this
                            print "Warning: conflicting strandedness information between BAM and GTF files. Read skipped." > "/dev/stderr"
                            print $0 > "/dev/stderr"
                        }else{
                            stats["reads_without_mismatches_filtered_final"]++
                            if (strand == "+" || strand2 == "forward"){
                                if (NODUPS_ARG == 1){
                                    if (and($2, 0x400)){
                                        print gene_id, b[1], 0, 0, 0, 0, b[7], 0, 0, 0, 0, b[13], 0, 0, 0, 0, b[19], "T" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                    }else{
                                        print gene_id, b[1], 0, 0, 0, 0, b[7], 0, 0, 0, 0, b[13], 0, 0, 0, 0, b[19], "F" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                    }
                                }else{
                                    print gene_id, b[1], 0, 0, 0, 0, b[7], 0, 0, 0, 0, b[13], 0, 0, 0, 0, b[19] > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                }
                            }else if (strand == "-" || strand2 == "reverse"){
                                if (NODUPS_ARG == 1){
                                    if (and($2, 0x400)){
                                        print gene_id, b[19], 0, 0, 0, 0, b[13], 0, 0, 0, 0, b[7], 0, 0, 0, 0, b[1], "T" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                    }else{
                                        print gene_id, b[19], 0, 0, 0, 0, b[13], 0, 0, 0, 0, b[7], 0, 0, 0, 0, b[1], "F" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                    }
                                }else{
                                    print gene_id, b[19], 0, 0, 0, 0, b[13], 0, 0, 0, 0, b[7], 0, 0, 0, 0, b[1] > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                }
                            }else{
                                print "Error: Cannot determine strandedness. Script terminated." > "/dev/stderr"
                                print $0 > "/dev/stderr"
                                exit 1 # remove later, this likely only problematic is slamdunk mode, affects only reads which were paired, but did not have a mate (pair orientation is used to determine strandedness).
                            }
                        }
                    }
                }
            }
        }else{
            stats["reads_with_mismatches"]++
            ec=cigar_extender($6)
            MM_extractor($10, $11, ec)
            md_ext=MD_extender(md)
            TCRA_producer(mm_read, mm_quali, md_ext)
            # if (MAPPER_ARG == "bbmap"){
            #     ed_value=xi
            # }
            if ($5 >= MAPQ_FILTER && ed_value <= NM_FILTER && xi >= MISMATCH_FILTER){
                stats["reads_with_mismatches_filtered"]++
                if (MODE_ARG == "slamdunk"){ # MP tags are made and added extra
                    MP_producer(md, ec, $10)
                    print $0,"XI:f:"xi,"TC:i:"tc,"RA:Z:"ra,mp
                }else if (MODE_ARG == "normal"){
                    print $0,"XI:f:"xi,"TC:i:"tc,"RA:Z:"ra
                    if (xf !~ /__ambiguous|__no_feature|_too_low_aQual|__not_aligned|__alignment_not_unique/){ # move this filter before tha main loop to avoid processing of reads that will be excluded anyways by HTseq
                        split(xf_value,a,"__") # only if annotation is appended with strandedness
                        gene_id = a[1]; strand = a[2]
                        split(ra,b,",")
                        if ((strand == "+" && strand2 == "reverse") || (strand == "-" && strand2 == "forward")){
                            print "Warning: conflicting strandedness information between BAM and GTF files. Read skipped." > "/dev/stderr"
                            print $0 > "/dev/stderr"
                        }else{
                            stats["reads_with_mismatches_filtered_final"]++
                            if (strand == "+" || strand2 == "forward"){
                                if (NODUPS_ARG == 1){
                                    if (and($2, 0x400)){
                                        print gene_id, b[1], b[2], b[3], b[4], b[6], b[7], b[8], b[9], b[11], b[12], b[13], b[14], b[16], b[17], b[18], b[19], "T" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                    }else{
                                        print gene_id, b[1], b[2], b[3], b[4], b[6], b[7], b[8], b[9], b[11], b[12], b[13], b[14], b[16], b[17], b[18], b[19], "F" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                    }
                                }else{
                                    print gene_id, b[1], b[2], b[3], b[4], b[6], b[7], b[8], b[9], b[11], b[12], b[13], b[14], b[16], b[17], b[18], b[19] > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                }
                            }else if (strand == "-" || strand2 == "reverse"){
                                if (NODUPS_ARG == 1){
                                    if (and($2, 0x400)){
                                        print gene_id, b[19], b[18], b[17], b[16], b[14], b[13], b[12], b[11], b[9], b[8], b[7], b[6], b[4], b[3], b[2], b[1], "T" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                        }else{
                                            print gene_id, b[19], b[18], b[17], b[16], b[14], b[13], b[12], b[11], b[9], b[8], b[7], b[6], b[4], b[3], b[2], b[1], "F" > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                        }
                                    }else{
                                        print gene_id, b[19], b[18], b[17], b[16], b[14], b[13], b[12], b[11], b[9], b[8], b[7], b[6], b[4], b[3], b[2], b[1] > OUTDIR_ARG "/" BAM_NAME "_slam_counts.txt"
                                    }
                                }else{
                                    print "Error: Cannot determine strandedness. Script terminated." > "/dev/stderr"
                                    print $0 > "/dev/stderr"
                                    exit 1 # remove later, this likely only problematic is slamdunk mode, affects only reads which were paired, but did not have a mate (pair orientation is used to determine strandedness).
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

END{
    finish_time = systime()
    print "Reads with matches or mismatches:", stats["reads_with_M_in_cig"] >> OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log"
    print "Reads with mismatches:", stats["reads_with_mismatches"] >> OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log"
    print "Reads with mismatches passing filters:", stats["reads_with_mismatches_filtered"] >> OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log"
    if (MODE_ARG == "normal"){
        print "Reads with mismatches passing filters (unique assignment):", stats["reads_with_mismatches_filtered_final"] >> OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log"
    }
    print "Reads without mismatches:", stats["reads_without_mismatches"] >> OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log"
    print "Reads without mismatches passing filters:", stats["reads_without_mismatches_filtered"] >> OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log"
    if (MODE_ARG == "normal"){
        print "Reads without mismatches passing filters (unique assignment):", stats["reads_without_mismatches_filtered_final"] >> OUTDIR_ARG "/" BAM_NAME "_rebound_stats.log"
    }
    print "S(L)AM converter ended @ " strftime() > "/dev/stderr"
    print (finish_time - start_time)/3600 , "Execution time in hours" > "/dev/stderr"
    print (finish_time - start_time)/60 , "Execution time in minutes" > "/dev/stderr"
}
