BEGIN {
    print "S(L)AM converter started @ " strftime() > "/dev/stderr"
    OFS="\t"
    FS="\t"
}

# --- FUNCTIONS --- #
function tag_finder(tag,    i){
    if (substr($0,1,1)!="@"){
        for (i=1; i<=NF+1; i++){
            if (i==NF+1) {
                print "___________________________________________________"
                print "Fatal error: unable to find " tag " at record " $0 # > "/dev/stderr"
                exit 1
            }else if (substr($i,1,5)!=tag){
                continue
            }else{
                return $i
            }
        }
    }
}

function MD_extender(md,  out,a,b) {
    out=""
    # print md, "md in fucntion"
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
    return out ####
}
# returns matched+mismatched portion (--> out) of the original read (read) as an array (--> cr)
function MM_extractor(read, ec,     a,b,j,i,k,out) {
    out=""
    split(read, b, "")
    split(ec, a, "")
    j=1
    for (i in a) {
        # print a[i], i
        if (a[i] ~ /[I]/) {
            j++
            continue
        }else if (a[i] ~ /[S]/) {
            j++
            continue
        }else if (a[i] ~ /[DN]/) {
            continue
        }else if (a[i] ~ /[M]/) {
            out = out b[j]
            j++
            continue
        }else{
            print "Fatal error: Unexpected character(s) in CIGAR string of read" $1
            exit 1
        }
    }
    return out
}
function cigar_extender(cygaro,   out, b,i,j){
    # expands CIGAR string into a longer string (out) with a single letter per alignment position; works for CIGARs with [SIDNM] operations
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
    # print "a", a
    n=patsplit(a, b, /[[:alpha:]^]/, seps)
    for(i=1; i<=n; i++){
        #print "sepsi",seps[i], i
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
        # print "sum_md",sum_md
        if (b[i]!="^" && b[i]!="h"){
            ref_mis[sum_md]=b[i]
        }
    }
    # for (i in ref_mis){
    #     print "ref_misi",ref_mis[i], i
    # }
    # print ec
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
            # print "j",j
            # print offset[i], i,"offset"
                # offset[i] = j
        }

        offset[i] = j

        # print i, j, i+j
        # print substr(read,i+j,1), "substr"
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
        # remove later, this should never be required
        # gsub(",,*" , "," , mp)
        # gsub("MP:Z:," , "MP:Z:" , mp)
        if (length(ref_mis) == 0){
            print "This should not have happened" > "/dev/stderr"
            exit 1
        }
    }
#     for (i in offset) {
#     print offset[i],i, "off"
# }
}

function TCRA_producer(MM_extracted, md_extended,    i){
    delete count; delete mmms
    split(MM_extracted, mm_split, "")
    split(md_extended, md_split, "")
    for (i=1; i<=length(mm_split); i++) {
        if (md_split[i] == ".") {
            mmms[i] = mm_split[i] mm_split[i]
        }else{
            mmms[i] = md_split[i] mm_split[i]
        }
    }
    for (i in mmms) {
        if       (mmms[i] == "AA"){
            count["AA"]++
        }else if (mmms[i] == "AC"){
            count["AC"]++
        }else if (mmms[i] == "AG"){
            count["AG"]++
        }else if (mmms[i] == "AT"){
            count["AT"]++
        }else if (mmms[i] == "AN"){
            count["AN"]++
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
        }else if (mmms[i] == "TA"){
            count["TA"]++
        }else if (mmms[i] == "TC"){
            count["TC"]++
        }else if (mmms[i] == "TG"){
            count["TG"]++
        }else if (mmms[i] == "TT"){
            count["TT"]++
        }else if (mmms[i] == "TN"){
            count["TN"]++
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
        }
    }
    # reverse strand + first in pair or forward strand + second in pair
    if ((and($2, 0x20) && and($2, 0x80)) || (and($2, 0x10) && and($2, 0x40))){
        tc = "TC:i:" count["TC"]+0
    #reverse strand + second in pair or forward strand + first in pair
    }else if ((and($2, 0x10) && and($2, 0x80)) || (and($2, 0x20) && and($2, 0x40))){
        tc = "TC:i:" count["AG"]+0
    }
    sp = ","
    ra = "RA:Z:" count["AA"]+0 sp count["AC"]+0 sp count["AG"]+0 sp count["AT"]+0 sp count["AN"]+0 sp count["CA"]+0 sp count["CC"]+0 sp count["CG"]+0 sp count["CT"]+0 sp count["CN"]+0 sp count["GA"]+0 sp count["GC"]+0 sp count["GG"]+0 sp count["GT"]+0 sp count["GN"]+0 sp count["TA"]+0 sp count["TC"]+0 sp count["TG"]+0 sp count["TT"]+0 sp count["TN"]+0 sp count["NA"]+0 sp count["NC"]+0 sp count["NG"]+0 sp count["NT"]+0 sp count["NN"]+0 sp

    if (count["AC"] + count["AG"] + count["AT"] + count["AN"] + count["CA"] + count["CG"] + count["CT"] + count["CN"] + count["GA"] + count["GC"] + count["GT"] + count["GN"] + count["TA"] + count["TC"] + count["TG"] + count["TN"] + count["NA"] + count["NC"] + count["NG"] + count["NT"] + count["NN"] > 0){
        xi = "XI:f:"(count["AA"] + count["CC"] + count["GG"] + count["TT"] + 0) / (count["AA"] + count["AC"] + count["AG"] + count["AT"] + count["AN"] + count["CA"] + count["CC"] + count["CG"] + count["CT"] + count["CN"] + count["GA"] + count["GC"] + count["GG"] + count["GT"] + count["GN"] + count["TA"] + count["TC"] + count["TG"] + count["TT"] + count["TN"] + count["NA"] + count["NC"] + count["NG"] + count["NT"] + count["NN"])
    }else{
        xi = "XI:f:0"
    }
}

function TCRA_producer_simple(read, ec,      a, count) {
    # use "simple" function also if reads were entirely clipped (cases of 100% overlap) == only S in CIGARS; this reads may actually be tossed completely, but they can be used for DE analysis when one wants to use exactly the same pool of reads for both, DE and SLAM-seq; in such sase, the reads would need to get the identi score of 1, otherwise will be filtered out
    # no matches and mismatches == full clip
    if ($6 ~ /^[^M]+$/){
        tc_simple = "TC:i:0"
        ra_simple = "RA:Z:0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
        xi_simple = "XI:f:0" # may be set to 1 if filtering of reads with no Ms in CIGARs is not desired
    }else{
            # CIGAR case with M only
        if ($6 ~ /^[^SDNI]+$/) {
            a = $10
            # CIGAR cases with I,N,D,S
        }else{
            a = MM_extractor($10, ec)
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
        tc_simple = "TC:i:0"
        sp = ","
        ra_simple = "RA:Z:" count["A"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["C"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["G"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["T"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["N"]+0 sp
        xi_simple = "XI:f:1"
    }
}

# ---- main routine ---- #

# get MD tag field if present
{
    if ($0 ~ /^.*MD:Z:.*$/){
    md = tag_finder("MD:Z:")
}
    # skip header and modify @RG by adding DS. SM can be supplied as argument. To uses files normally produced by nf-core pipeline. Alternatively, STAR + Picard MarkDuplicates
    if ($1 ~ /^@/){
        header=1
        if ($1 ~ /^@RG/){
            #rg_line=$1
        #    for (i=2; i<=NF; i++){
            #    rg_line=rg_line","$i
            #}
            suf="DS:{'sequenced':0,'mapped':0,'filtered':0,'mqfiltered':0,'idfiltered':0,'nmfiltered':0,'multimapper':0,'dedup':0,'snps':0,'annotation':'','annotationmd5':''}"
            #rg_line=$1 ", ID:" BAM_NAME ", SM:" BAM_NAME ":NA:NA," suf
            print $1,"ID:"BAM_NAME,"SM:"BAM_NAME":NA:NA", suf
        }else{
               print $0
        }
    }else if ($1 !~ /^@/ && header==1){
        header=0
        print "@PG","ID:rebound.sh","VN:1.0","CL:"REBOUND_COMMAND
        if (and($2, 0x4)){
            print $0, "XI:f:0"
        }else if (gensub(/MD:Z:/, "", "g", md) ~ /^[0-9]*(\^[ACTGN]+[0-9]+)*$/ || $6 ~ /^[^M]+$/){
            ec=cigar_extender($6)
            TCRA_producer_simple($10, ec)
            print $0,xi_simple,tc_simple,ra_simple
        }else{
            ec=cigar_extender($6)
            MM_extracted=MM_extractor($10, ec)
            md_ext=MD_extender(md)
            MP_producer(md, ec, $10)
            TCRA_producer(MM_extracted, md_ext)
            print $0,xi,tc,ra,mp
        }
    }else if (and($2, 0x4)){
        print $0, "XI:f:0"
    # }else if ((and($2, 0x10) && and($2, 0x80)) || (and($2, 0x20) && and($2, 0x40))){
    #     print $0, "XI:f:0"
    }else if (gensub(/MD:Z:/, "", "g", md) ~ /^[0-9]*(\^[ACTGN]+[0-9]+)*$/ || $6 ~ /^[^M]+$/){
        ec=cigar_extender($6)
        TCRA_producer_simple($10, ec)
        print $0,xi_simple,tc_simple,ra_simple
    }else{
        ec=cigar_extender($6)
        MM_extracted=MM_extractor($10, ec)
        md_ext=MD_extender(md)
        MP_producer(md, ec, $10)
        TCRA_producer(MM_extracted, md_ext)
        print $0,xi,tc,ra,mp
        }
    }

END{
    print "S(L)AM converter ended @ " strftime() > "/dev/stderr"
}
